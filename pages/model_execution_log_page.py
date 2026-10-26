"""Model Execution Log page.

Shown while a model run is in progress. The whole page is one animated card,
`widgets.omics_loader.OmicsLoader`: the run's input files, an illustrative
graph that builds up as the run proceeds, a console streaming the backend's
real log lines, and the current stage with a five-part progress bar.

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) with "Model Logs"
      marked as the active nav item.
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Model Running" marked as the active workflow tab.
    - The loader card, which this page only feeds: log lines, progress, and
      whether the run completed or failed.

Interactions with other pages:
    - `on_upload_clicked` navigates back to `DatasetInitializationPage`.
    - `on_model_visualization_clicked` navigates to `ModelVisualizationPage`.
    - `on_run_complete` is invoked once the watched run has finished
      successfully and the loader has played through to "Done", and should
      navigate to the results page.
    - `on_results_clicked` navigates back to `FinalResultsPage` from the
      header's "Results" tab.
    All four callbacks are supplied and wired by `UILauncher.py`.

Call `start_watching(run_info)` (with the dict returned by
`client.api_client.start_run`, merged with the upload response) after
constructing this page to begin polling `client.workers.RunStatusPoller`
for live log lines and progress.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFrame, QHBoxLayout, QMessageBox, QScrollArea, QSizePolicy, QVBoxLayout, QWidget

from client.workers import RunStatusPoller
from styles.theme import SECTION_SPACING
from widgets.navigation import (
    build_header_bar,
    build_sidebar,
    make_primary_cta_button,
    make_sidebar_nav_button,
    make_top_tab,
)
from widgets.omics_loader import OmicsLoader

# What each of the loader's file cards shows for the final-model backend. The
# graph supplies the drug-target and driver-mutation edges the pair features
# are computed from; each omics layer is its own PCA file. The drug-response
# card shows the uploaded file's own name (see `start_watching`).
_RUN_FILES = {
    "graph": "hetero_graph.pt",
    "expression": "transcriptomics_pca.csv",
    "mutation": "genomics_pca.csv",
    "proteomics": "proteomics_pca.csv",
    "model": "final model · 5 seeds",
}


class ModelExecutionLogPage(QWidget):
    """Live view of an in-progress model run."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_upload_clicked: Callable[[], None] | None = None,
        on_model_visualization_clicked: Callable[[], None] | None = None,
        on_run_complete: Callable[[], None] | None = None,
        on_results_clicked: Callable[[], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_upload_clicked: Invoked when the header's "Upload" tab is
                clicked; should navigate back to the dataset upload page.
            on_model_visualization_clicked: Invoked when the sidebar's
                "Model Visualization" link is clicked; should navigate to
                the model visualization page.
            on_run_complete: Invoked once the watched run finishes
                successfully; should navigate to the results page.
            on_results_clicked: Invoked when the header's "Results" tab is
                clicked; should return to the results page of the run
                already shown there, without reloading it.
        """
        super().__init__(parent)
        self.setObjectName("ModelExecutionLogPage")
        self._on_upload_clicked = on_upload_clicked
        self._on_model_visualization_clicked = on_model_visualization_clicked
        self._on_run_complete = on_run_complete
        self._on_results_clicked = on_results_clicked

        self._poller: RunStatusPoller | None = None
        self._current_run_id: str | None = None
        self._loader = OmicsLoader()
        self._loader.finished.connect(self._on_loader_finished)

        self._build_ui()

    # -- live updates ---------------------------------------------------

    def start_watching(self, run_info: dict) -> None:
        """Begin polling and displaying live status for a newly started run.

        Args:
            run_info: The merged upload+run-start response from
                `client.workers.UploadWorker` (must include `run_id`;
                `filename` and `target_cell_line` are shown when present).
        """
        if self._poller is not None:
            self._poller.stop()
            try:
                self._poller.statusUpdate.disconnect(self._on_status_update)
            except (RuntimeError, TypeError):
                pass

        self._current_run_id = run_info["run_id"]
        files = {**_RUN_FILES, "response": run_info.get("filename") or "uploaded dataset"}
        self._loader.start(files, str(run_info.get("target_cell_line") or ""))

        self._poller = RunStatusPoller(run_info["run_id"], parent=self)
        self._poller.statusUpdate.connect(self._on_status_update)
        self._poller.start()

    def _on_status_update(self, payload: dict) -> None:
        # A poller for a run we've since replaced (e.g. a new upload started
        # while this one was still finishing up) may still have one update
        # in flight when it's told to stop; ignore anything not for the run
        # currently being watched so a late "completed" from an old run
        # can't jump the page to results out from under the user.
        if payload.get("run_id") != self._current_run_id:
            return

        for line in payload.get("new_log_lines", []):
            self._loader.log(line)
        self._loader.set_progress(payload.get("progress_percent", 0.0) / 100.0)

        status = payload.get("status", "running")
        if status == "completed":
            # The loader plays out its remaining steps, then emits `finished`.
            self._loader.complete()
        elif status == "failed":
            message = payload.get("error_message") or "Unknown error."
            self._loader.fail(message)
            QMessageBox.critical(self, "Model Run Failed", message)

    def _on_loader_finished(self) -> None:
        if self._on_run_complete is not None:
            self._on_run_complete()

    def _build_ui(self) -> None:
        """Lay out the sidebar, header, and the loader card."""
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        root.addWidget(self._build_sidebar(), 0)

        content_shell = QFrame()
        content_shell.setObjectName("RootShell")
        content_shell_layout = QVBoxLayout(content_shell)
        content_shell_layout.setContentsMargins(0, 0, 0, 0)
        content_shell_layout.setSpacing(0)
        content_shell.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        content_shell_layout.addWidget(self._build_header(), 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        # Tighter than PAGE_MARGIN so the whole card fits the default window without scrolling.
        content_layout.setContentsMargins(SECTION_SPACING, SECTION_SPACING, SECTION_SPACING, SECTION_SPACING)
        content_layout.addWidget(self._loader, 1)

        scroll.setWidget(content)
        content_shell_layout.addWidget(scroll, 1)
        root.addWidget(content_shell, 1)

    def _build_sidebar(self) -> QFrame:
        """Build the shared sidebar with "Model Logs" as the active nav item."""
        nav_widgets = [
            # !Hidden: ModelVisualizationPage is a decorative placeholder
            # (fake canvas + hardcoded run stats), hidden from navigation
            # until it's wired to real data.
            # make_sidebar_nav_button(
            #     "Model Visualization",
            #     "insights",
            #     callback=self._on_model_visualization_clicked,
            # ),
            make_sidebar_nav_button("Model Logs", "terminal", active=True),
        ]
        footer_widgets = [make_sidebar_nav_button("Support", "help_outline")]
        return build_sidebar(
            nav_widgets,
            footer_widgets,
            cta_widget=make_primary_cta_button(),
        )

    def _build_header(self) -> QFrame:
        """Build the shared header bar with "Model Running" as the active tab."""
        tabs = [
            make_top_tab("Upload", callback=self._on_upload_clicked),
            make_top_tab("Model Running", active=True),
            make_top_tab("Results", callback=self._on_results_clicked),
        ]
        return build_header_bar(tabs)
