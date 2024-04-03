"""
The controller portion of the PypeIt Setup GUI.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import traceback
import signal
import sys
from datetime import datetime
import re
import io
from pathlib import Path
from functools import partial
from contextlib import contextmanager

# TODO: datetime.UTC is not defined in python 3.10.  Remove this when we decide
# to no longer support it.
try:
    __UTC__ = datetime.UTC
except AttributeError as e:
    from datetime import timezone
    __UTC__ = timezone.utc

from qtpy.QtCore import QCoreApplication, Signal, QMutex
from qtpy.QtCore import QObject, Qt, QThread
from qtpy.QtGui import QKeySequence
from qtpy.QtWidgets import QAction

from pypeit.setup_gui.view import SetupGUIMainWindow, DialogResponses
from pypeit.setup_gui.text_viewer import TextViewerWindow
from pypeit.setup_gui.dialog_helpers import prompt_to_save, display_error, FileDialog, FileType
from pypeit.setup_gui.model import PypeItSetupGUIModel, ModelState
from pypeit import msgs

from pypeit.display import display
from pypeit import io as pypeit_io


# For using Qt Mutexes and "with"
@contextmanager
def lock_qt_mutex(mutex):
    mutex.lock()
    try:
        yield mutex
    finally:
        mutex.unlock()

class OpCanceledError(Exception):
    """Exception thrown when a background operation has been canceled."""
    def __init__(self):
        super().__init__()

class OperationThread(QThread):
    """Thread to run a background operation."""

    completed = Signal(bool, tuple)
    """Signal send the operation has completed."""


    def __init__(self):
        super().__init__()
        self._operation = None
        self._max_progress = None
        self._mutex = QMutex()

    def run(self):
        """Runs an operation in a background thread."""
        canceled = False
        exc_info = (None, None, None)
        try:
            self._operation.run()            
        except OpCanceledError:
            canceled=True
        except Exception:
            exc_info=sys.exc_info()

        self.completed.emit(canceled, exc_info)        

    def _cancel_op(self):
        """Cancels an in progress background operation when the user cancels the progress dialog."""
        self.requestInterruption()

    def _op_progress(self, max_progress, progress_message=None):
        """Signal handler that receives progress information from a background operation
        as it proceeds. It passes this to the view to increase the value showing in the
        progress dialog."""

        with lock_qt_mutex(self._mutex):
            create_progress = False
            mp = None
            if self._operation is not None:
                if self._max_progress is None:
                    if max_progress is not None:
                        self._max_progress = max_progress
                        create_progress = True
                mp = self._max_progress
        
        if create_progress:
            # If we've just initialized the max progress, create the progress dialog
            SetupGUIController.main_window.create_progress_dialog(self._operation.name, max_progress, self._cancel_op)
            
        # Ignore the progress if there's no max progress yet
        if mp is not None:            
            SetupGUIController.main_window.show_operation_progress(increase=1, message = progress_message)

    def _op_complete(self, canceled, exc_info):
        """Signal handler that is notified when a background operation completes.
        
        Args:
            canceled (bool): Whether or not the operation was canceled.
            exc_info (tuple): The exception information if the operation failed. None if it succeeded
        """
        msgs.info("Op complete")
        with lock_qt_mutex(self._mutex):
            if self._operation is not None:
                operation = self._operation
                self._operation = None
                self._max_progress = None
            else:
                operation = None

        if operation is not None:
            operation.progressMade.disconnect(self._op_progress)
            SetupGUIController.main_window.operation_complete()
            operation.postRun(canceled, exc_info)            
    

    def startOperation(self, operation):
        """
        Start a background operation.
        
        Args:
            operation (MetadataOperation): The MetadataOperation to start in the background thread.
        """
        self._operation = operation
        if operation.preRun():
            operation.progressMade.connect(self._op_progress, type=Qt.QueuedConnection)
            self.completed.connect(self._op_complete, type=Qt.QueuedConnection)
            self.start()

class MetadataOperation(QObject):

    """Base class for Metadata operations that take long enough that they should take place in a background thread.
    
    Args:
        model (PypeItSetupGUIModel): The PypeIt Setup GUI's model object.
    """

    progressMade = Signal(int, str)
    """Signal emitted emit when progress has been made. This will be reflected in the view's progress dialog."""


    def __init__(self, name, model):
        super().__init__()
        self._model=model
        self.name = name
        self._max_progress = None

    def preRun(self):
        """
        Perform setup required before running the operation. This involves watching the log
        for files being added to the metadata.
        """
        building_metadata_re = re.compile("Building metadata for (\d+) ")
        self._model.log_buffer.watch("building_metadata", building_metadata_re, self._buildingMetadata)
        
        added_metadata_re = re.compile("Adding metadata for (.*)$")
        self._model.log_buffer.watch("added_metadata", added_metadata_re, self._addedMetadata)
        self._model.closeAllFiles()
        return True

    def _buildingMetadata(self, name, match):
        """Callback used to find the total number of files being read when building metadata."""
        self._max_progress = int(match.group(1))
        msgs.info(f"Found max progress {self._max_progress}")

    def _addedMetadata(self, name, match):
        """Callback used to report progress on reading files when building metadata."""
        if QThread.currentThread().isInterruptionRequested():
            raise OpCanceledError()

        self.progressMade.emit(self._max_progress, match.group(1))


    def postRun(self, canceled, exc_info):
        """Clean up steps after the operations has run.
        
        Args:
            canceled (bool):  True if the operation was canceled. 
            exc_info (tuple): The exception information (as returned by sys.exc_info()) for any errors that occurred.
        """
        self._model.log_buffer.unwatch("added_metadata")
        self._model.log_buffer.unwatch("building_metadata")

        if exc_info[0] is not None:
            traceback_string = "".join(traceback.format_exception(*exc_info))
            msgs.warn(f"Failed to {self.name.lower()}:\n" + traceback_string)
            display_error(SetupGUIController.main_window, f"Failed to {self.name.lower()} {exc_info[0]}: {exc_info[1]}")
            self._model.reset()
        elif canceled:
            self._model.reset()

    def run(self):
        """Performs the steps of the MetadataOperation. This is an abstract method overridden by child classes."""
        pass

class SetupOperation(MetadataOperation):
    """ Background operation to run pypeit setup.

        Args:
            model (PypeItSetupGUIModel): The PypeIt Setup GUI's model object.
    """
    def __init__(self, model):
        super().__init__("Run Setup", model)

    def run(self):
        """
        Runs pypeit_setup on the current raw data directories.
        """
        self._model.run_setup()

class OpenFileOperation(MetadataOperation):
    """
    Background operation to open a PypeIt file

    Args:
        model (PypeItSetupGUIModel):
            The PypeIt Setup GUI's model object.
        file (): 
            The file to open.
    """

    def __init__(self, model, file):
        super().__init__("Open File", model)
        self._file = file

    def run(self):
        """
        Opens a pypeit file and reads metadata for all of the files in it.
        """
        self._model.open_pypeit_file(self._file)

class MetadataReadOnlyAction(QAction):
    """An action (caused by the right click menu in the GUI, a button, or keyboard short cut) that is read only
    and therefore can be performed on the ObsLog.
    
    Arguments:
        controller(PypeItMetadataController): 
            The controller for the PypeItMetadataModel/PypeItMetadataView MVC triplet.
        menu_text(str):
            The text name for the menu/button that triggers the action.
        handler (Callable):
            The signal handler to enact the action. This receives the "triggered" event
            from the parent class.
        shortcut (PySide2.QtGui.QKeySequence.StandardKey, Optional):
            The keyboard shortcut to initiate the action.
    """
    def __init__(self, controller, menu_text, handler, shortcut=None):
        super().__init__(menu_text)
        self.triggered.connect(handler)
        if shortcut is not None:
            self.setShortcut(shortcut)
        self._controller=controller

    def updateEnabledStatus(self):
        """Enable/disable the action based on whether any metadata rows are selected."""
        if self._controller._view is not None and len(self._controller._view.selectedRows()) > 0:
            self.setEnabled(True)
        else:
            self.setEnabled(False)

class MetadataWriteAction(QAction):
    """An action (caused by the right click menu in the GUI, a button, or keyboard short cut) that can change
    the file metadata and therefore can only be performed on a PypeItFileModel.
    
    Arguments:
        controller(PypeItMetadataController): 
            The controller for the PypeItMetadataModel/PypeItMetadataView MVC triplet.
        menu_text(str):
            The text name for the menu/button that triggers the action.
        handler (Callable):
            The signal handler to enact the action. This receives the "triggered" event
            from the parent class.
        shortcut (PySide2.QtGui.QKeySequence.StandardKey, Optional):
            The keyboard shortcut to initiate the action.
    """

    def __init__(self, controller, menu_text, handler, shortcut=None):
        super().__init__(menu_text)
        self.triggered.connect(handler)
        if shortcut is not None:
            self.setShortcut(shortcut)
        self._controller=controller

    def updateEnabledStatus(self):
        """Enable/disable the action based on whether any metadata rows are selected."""

        if self._controller._is_pypeit_file:
            if self._controller._view is not None and len(self._controller._view.selectedRows()) > 0:
                self.setEnabled(True)
            else:
                self.setEnabled(False)
        else:
            self.setEnabled(False)

class MetadataPasteAction(QAction):
    """An action (caused by the right click menu in the GUI, a button, or keyboard short cut) for pasting
    metadata into a PypeItMetadataModel object.
    
    Arguments:
        controller(PypeItMetadataController): 
            The controller for the PypeItMetadataModel/PypeItMetadataView MVC triplet.
        menu_text(str):
            The text name for the menu/button that triggers the action.
        handler (Callable):
            The signal handler to enact the action. This receives the "triggered" event
            from the parent class.
        shortcut (PySide2.QtGui.QKeySequence.StandardKey, Optional):
            The keyboard shortcut to initiate the action.
    """

    def __init__(self, controller, menu_text, handler, shortcut=None):
        super().__init__(menu_text)
        self.triggered.connect(handler)
        if shortcut is not None:
            self.setShortcut(shortcut)
        self._controller=controller

    def updateEnabledStatus(self):
        """Enable/disable the action based on whether any metadata rows are selected AND
        there are rows to paste in the clipboard."""

        if self._controller._is_pypeit_file:
            if SetupGUIController.model.clipboard.rowCount() > 0:
                self.setEnabled(True)
            else:
                self.setEnabled(False)
        else:
            self.setEnabled(False)

class PypeItMetadataController(QObject):
    """A Controller object for performing actions iniitiated by the user file metadata.
    Part of a MVC triplet involving PypeItMetadataModel/PypeItMetadataController/PypeItMetadataView.

    Args:
        model (PypeItMetatadataModel): 
            The model this controller acts with.

        is_pypeit_file (bool): 
            True if the model is for a PypeItFileModel (that is writeable model), False if it is 
            from a PypeItObsLog model (read only)         
    """
    def __init__(self, model, is_pypeit_file):
        super().__init__()
        self._model = model
        self._view = None
        self._is_pypeit_file = is_pypeit_file
        self._windows = {}
        self.next_window_id = 1

        # Build actions
        self._action_list = [MetadataReadOnlyAction(self, "View File",   self.view_file),
                             MetadataReadOnlyAction(self, "View Header", self.view_header),
                             MetadataReadOnlyAction(self, "Copy",        self.copy_metadata_rows,       shortcut=QKeySequence.StandardKey.Copy),
                                MetadataWriteAction(self, "Cut",         self.cut_metadata_rows,        shortcut=QKeySequence.StandardKey.Cut),
                                MetadataPasteAction(self, "Paste",       self.paste_metadata_rows,      shortcut=QKeySequence.StandardKey.Paste),
                                MetadataWriteAction(self, "Comment Out", self.comment_out_metadata_rows),
                                MetadataWriteAction(self, "Uncomment",   self.uncomment_metadata_rows),
                                MetadataWriteAction(self, "Delete",      self.delete_metadata_rows,     shortcut=QKeySequence.StandardKey.Delete) ]
        SetupGUIController.model.clipboard.modelReset.connect(self.updatedEnabledActions)
        self.updatedEnabledActions()

    def getActions(self, parent):
        """Returns the actions that this controller supports.
        
        Returns: (list of QAction): List of the actions that can be performed on the PypeItMetadataModel.
        """
        return self._action_list
            
    def setView(self, view):
        """Set the view that is responsible for displaying and receiving input for
        the PypeItMetadataModel.
        
        Args:
            view (PypeItMetadataView): The view.
        """
        self._view = view
        view.selectionUpdated.connect(self.updatedEnabledActions)

    def updatedEnabledActions(self):
        """Updates which actions are enabled/disabled."""
        for action in self._action_list:
            action.updateEnabledStatus()
                    
    def view_file(self):
        """View the selected files in the metadata using Ginga."""
        row_indices = self._view.selectedRows()
        if len(row_indices) > 0:

            # Make sure ginga is available
            try:
                display.connect_to_ginga(raise_err=True, allow_new=True)
            except Exception as e:
                display_error(SetupGUIController.main_window, f"Could not start ginga to view FITS files: {e}")
                msgs.warn(f"Failed to connect to ginga:\n" + traceback.format_exc())

            
            # Display each file in its own ginga tab
            for indx in row_indices:
                metadata_row = self._model.metadata[indx]
                file = Path(metadata_row['directory'], metadata_row['filename'])
                if not file.exists():
                    display_error(SetupGUIController.main_window, f"Could not find {file.name} in {file.parent}.")
                    return

                try:
                    img = self._model.spectrograph.get_rawimage(str(file), 1)[1]
                except Exception as e:
                    display_error(SetupGUIController.main_window, f"Failed to read image {file.name}: {e}")
                    msgs.warn(f"Failed get raw image:\n" + traceback.format_exc())

                try:
                    display.show_image(img, chname = f"{file.name}")
                except Exception as e:
                    display_error(SetupGUIController.main_window, f"Failed to send image {file.name} to ginga: {e}")
                    msgs.warn(f"Failed send image to ginga:\n" + traceback.format_exc())

    def view_header(self):
        """ Display the header of one or more selected files in the metadata.
        """
        row_indices = self._view.selectedRows()
        if len(row_indices) > 0:            
            # Display each file in its window
            for indx in row_indices:
                metadata_row = self._model.metadata[indx]
                file = Path(metadata_row['directory'], metadata_row['filename'])
                header_string_buffer = io.StringIO()
                try:
                    with pypeit_io.fits_open(file) as hdul:
                        for i, hdu in enumerate(hdul):
                            print(f"\n\n# HDU {i} Header from {file.name}\n",file=header_string_buffer)
                            hdu.header.totextfile(header_string_buffer)
                except Exception as e:
                    display_error(SetupGUIController.main_window, f"Failed to read header from file {file.name} in {file.parent}: {e}")
                    msgs.warn(f"Failed to read header from {file}:\n" + traceback.format_exc())
                    return
                header_string_buffer.seek(0)
                window = TextViewerWindow(title=f"{file.name} Header", width=80, height=50,start_at_top=True, filename=file.parent / (file.name+".txt"), text_stream=header_string_buffer)
                self._windows[self.next_window_id] = window
                window.closed.connect(partial(self._window_closed, id=self.next_window_id))
                self.next_window_id+=1
                window.show()

    def _window_closed(self, id : int) -> None:
        """Clean up when a header viewer window is closed"""
        del self._windows[id]

    def copy_metadata_rows(self):
        """Copy metadata rows into the clipboard."""
        if self._view is not None:
            row_indices = self._view.selectedRows()
            msgs.info(f"Copying {len(row_indices)} rows to the clipboard.")
            if len(row_indices) > 0:
                row_model = self._model.createCopyForRows(row_indices)
                SetupGUIController.model.clipboard = row_model
                return True
        else:
            msgs.warn("Copy from controller with no view")
        return False

    def cut_metadata_rows(self):
        """Move metadata rows from the PypeItMetadataModel to the clipboard."""
        if self.copy_metadata_rows():
            return self.delete_metadata_rows()
        return False


    def paste_metadata_rows(self):
        """Insert rows from the clipboard into the PypeItMetadataModel"""
        clipboard = SetupGUIController.model.clipboard
        if clipboard.rowCount() > 0:
            try:
                msgs.info(f"Pasting {clipboard.rowCount()} rows")
                self._model.pasteFrom(clipboard)
            except Exception as e:
                traceback_string = "".join(traceback.format_exc())
                msgs.warn(f"Failed to paste metadata rows:\n" + traceback_string)
                display_error(SetupGUIController.main_window, f"Could not paste rows to this PypeIt file: {e}")


    def comment_out_metadata_rows(self):
        """Comment out one or more selected metadata rows."""
        if self._view is not None:
            row_indices = self._view.selectedRows()
            msgs.info(f"Commenting out {len(row_indices)} rows.")
            if len(row_indices) > 0:
                self._model.commentMetadataRows(row_indices)
    
    def uncomment_metadata_rows(self):
        """Uncomment previously commented out selected metadata rows."""
        if self._view is not None:
            row_indices = self._view.selectedRows()
            msgs.info(f"Unommenting out {len(row_indices)} rows.")
            if len(row_indices) > 0:
                self._model.uncommentMetadataRows(row_indices)

    def delete_metadata_rows(self):
        """Remove one or more selected rows from the PypeItMetadataModel."""
        if self._view is not None:
            row_indices = self._view.selectedRows()
            msgs.info(f"Removing {len(row_indices)} rows.")
            if len(row_indices) > 0:               
                self._model.removeMetadataRows(row_indices)
                return True
        else:
            msgs.warn("Copy from controller with no view")
        return False



class PypeItObsLogController(QObject):
    """PypeItObsLogController responsible for responding to user interaction
    as part of a MVC triplet with PypeItObsLogView and PypeItObsLogModel

    Args:
        main_window (:obj:`UserPromptDelegate`): A view object that can prompt the user.
        model (:obj:`PypeItObsLogModel`): The model for the obs log.
        operation_thread (:obj:`pypeit.setup_gui.controller.SetupGUIController`): The main Setup GUI controller.
    """

    def __init__(self, model, setup_gui_controller):
        super().__init__()
        self._model = model
        self.setup_gui_controller = setup_gui_controller

    def setModel(self, model):
        self._model = model

    def getMetadataController(self, model):
        return PypeItMetadataController(model, is_pypeit_file=False)


    def setSpectrograph(self, spectrograph_name):
        self._model.set_spectrograph(spectrograph_name)

        if self._model.state != ModelState.NEW:
            # Re-run setup with the new spectrograph
            self.setup_gui_controller.run_setup()

        else:
            self._model.set_spectrograph(spectrograph_name)

    def removePaths(self, rows):
        # Remove paths in reverse order, so that indexes don't change when a row is removed
        for row in sorted(rows, reverse=True):
            self._model.paths_model.removeRow(row)

    def addNewPath(self, new_path):
        """Add a new path to the observation log"""
        msgs.info(f"Adding new path {new_path}")
        self._model.add_raw_data_directory(new_path)

class PypeItFileController(QObject):
    """PypeItFileController responsible for responding to user interaction
    as part of a MVC triplet with PypeItFileView and PypeItFileModel

    Args:
        main_window (:obj:`UserPromptDelegate`): A view object that can prompt the user.
        model (:obj:`PypeItFileModel`): The model for the obs log.
    """

    def __init__(self, model):
        self._model = model

    def setModel(self, model):
        self._model = model

    def getMetadataController(self, model):
        return PypeItMetadataController(model, is_pypeit_file=True)


class SetupGUIController(QObject):
    """Controller for the PypeIt setup gui. It is responsible for initializing the GUI,
    and performing actions requested by the user.
    
    Args:
        args (:class:`argparse.Namespace`): The non-QT command line arguments used to start the GUI.
    """


    main_window = None
    model = PypeItSetupGUIModel()

    def __init__(self, args):
        super().__init__()

        QCoreApplication.setOrganizationName("PypeIt")
        QCoreApplication.setApplicationName("SetupGUI")
        QCoreApplication.setOrganizationDomain("pypeit.readthedocs.io")

        if args.logfile is not None:
            logpath = Path(args.logfile)
            if logpath.exists():
                timestamp = datetime.now(__UTC__).strftime("%Y%m%d-%H%M%S")
                old_log=logpath.parent / (logpath.stem + f".{timestamp}" + logpath.suffix)
                logpath.rename(old_log)
                
        self.model.setup_logging(args.logfile, verbosity=args.verbosity)
        if args.spectrograph is not None:
            self.model.obslog_model.set_spectrograph(args.spectrograph)
        if args.root is not None:
            for root_dir in args.root:
                self.model.obslog_model.add_raw_data_directory(root_dir)
        if args.spectrograph is not None and args.root is not None:
            self.run_setup_at_startup = True
        else:
            self.run_setup_at_startup = False

        self.model.obslog_model.default_extension = args.extension
        SetupGUIController.main_window = SetupGUIMainWindow(self.model, self)
        self.operation_thread = OperationThread()


    def getObsLogController(self, model):
        """Create the PypeItObsLogController as part of a MVC triplet.

        Args:
            model (:obj:`PypeItObsLogModel`): The model for the obs log.
            main_window (:obj:`SetupGUIMainWindow`): The mainwindow of the setup GUI.
        """
        return PypeItObsLogController(model, self)
    
    def getPypeItFileController(self, model):
        """Create the PypeItObsLogController as part of a MVC triplet.

        Args:
            model (:obj:`PypeItFileModel`): The model for the obs log.
        """
        return PypeItFileController(model)

    def start(self, app):
        """
        Starts the PypeItSetupGUi event loop. Exits the GUI when the GUI is closed.

        Args:
            app (QApplication): The Qt application object for the GUI. The caller is expected
                                to pass any Qt specific command line arguments to this object
                                before calling start(). 
        """
        self.main_window.show()
        if self.run_setup_at_startup:
            self.run_setup()

        # QT Gobbles up the Python Ctrl+C handling, so the default PypeIt Ctrl+C handler won't
        # be called. So we reset it to the OS default
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        sys.exit(app.exec_())

    def save_all(self):
        """"
        Save all unique configurations as pypeit files. Called in response to the user
        clicking the "Save All" button.
        """
        try:
            response = DialogResponses.CANCEL
            location = None
            for file_model in self.model.pypeit_files.values():
                if file_model.save_location is None:
                    if location is None or response != DialogResponses.ACCEPT_FOR_ALL:
                        dialog = FileDialog.create_save_location_dialog(self.main_window, file_model.name_stem, True)
                        response = dialog.show()
                        if response == DialogResponses.CANCEL:
                            # Cancel was pressed
                            return
                        location = dialog.selected_path

                file_model.save_location = location
                file_model.save()

        except Exception as e:
            display_error(self.main_window, str(e))

    def save_one(self):
        """ Saves the currently selected configuration as a pypeit file. Called in response to the user
        clicking the "Save Tab" button."""
        try:
            config_name = self.main_window.tab_widget.currentWidget().name
            self._save_tab(config_name)
        except Exception as e:
            display_error(self.main_window, str(e))
    
    def _save_tab(self, config_name):
        msgs.info(f"Saving config {config_name}")
        file_model = self.model.pypeit_files[config_name]
        if file_model.save_location is None:
            dialog = FileDialog.create_save_location_dialog(self.main_window, config_name)
            if dialog.show() == DialogResponses.CANCEL:
                return
            else:
                file_model.save_location = dialog.selected_path
        file_model.save()

    def clear(self):
        """Resets the GUI to it's initial state. Called in response to the user
        clicking the "Clear" button. This will prompt the user if they want to save
        any unsaved changes first."""

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.model.reset()

    def createNewFile(self, source_file_name, selectedRows):
        try:
            self.model.createNewPypeItFile(source_file_name, selectedRows)
        except Exception as e:
            display_error(self.main_window, f"Failed to create new tab {e.__class__.__name__}: {e}")
            msgs.warn(f"Failed to create new tab.")
            msgs.warn(traceback.format_exc())

    def close(self, file_model):

        try:
            if file_model.state == ModelState.CHANGED:
                response = prompt_to_save(self.main_window)
                if response == DialogResponses.SAVE:
                    self._save_tab(file_model.name_stem)
                elif response == DialogResponses.CANCEL:
                    return False
        except Exception as e:
            display_error(self.main_window, str(e))
            return False

        self.model.removeFile(file_model.name_stem)
        return True


    def exit(self):
        """Exits the GUI. Called in response to the user clicking the "Exit" button.
        This will prompt the user if they want to save any unsaved changes first."""

        # Prompt to save unsaved changes
        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        sys.exit(0)

    def run_setup(self):
        """Runs setup on the currently selected raw data directories. 
           Called in response to the user the clicking the "Run Setup" or changing the spectrograph.
           This will prompt the user if they want to save any unsaved changes first and
           start the operation in a background thread.
        """

        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        self.operation_thread.startOperation(SetupOperation(self.model))

    def createNewPypeItFile(self):
        # First figure out the name
        # TODO, use New File convention as above, which forces "save as" when saving a new file,
        # or current default single letter naming scheme?
        # If we use letters, use base 26 numbers to support more than 1 letter?
        #new_name = "New File"
        #i = 1
        #while new_name in self.pypeit_files.keys():
        #    i += 1
        #    new_name = f"New File{i}"


        if len(self.model.pypeit_files) == 0:
            # No configs, just add "A"
            new_name = "A"
        else:
            largest_name = max(self.model.pypeit_files.keys())
            if largest_name == 'z':
                raise ValueError("Failed to create new setup because there are too many loaded.")
            new_name = chr(ord(largest_name)+1)

        self.model.createEmptyPypeItFile(new_name)


    def open_pypeit_file(self):
        """Opens a PypeIt file. Called in response to the user the clicking the "Open" button. 
           This method prompts the user to discard or save any current changes,
           prompts the user for a pypeit to open, and opens it.
        """
        # Prompt to save current changes if needed
        if self.model.state == ModelState.CHANGED:
            response = prompt_to_save(self.main_window)
            if response == DialogResponses.SAVE:
                self.save_all()
            elif response == DialogResponses.CANCEL:
                return

        open_dialog = FileDialog.create_open_file_dialog(self.main_window, "Select PypeIt File", file_type=FileType("PypeIt input files",".pypeit"))
        result = open_dialog.show()
        if result != DialogResponses.CANCEL:
            self.operation_thread.startOperation(OpenFileOperation(self.model, open_dialog.selected_path))
