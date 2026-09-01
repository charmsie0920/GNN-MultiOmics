"""Dataset Initialization page.

The first step in the workflow: lets a user stage a dataset upload, verify a
target drug, select required cancer/omics types, and configure an analysis
profile before kicking off a model run.

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`).
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Upload" marked as the active tab and a "Results" dropdown menu
      (`QToolButton`) standing in for a future results-history picker.
    - An upload dropzone card and a "Target Sample" card, the two functional
      pieces of this page — everything else (target drug verification,
      cancer/omics-type selectors, analysis profile) was cosmetic-only and
      has been removed.

Interactions with other pages:
    - `on_initialize_upload` navigates to `ModelExecutionLogPage`, wired by
      `main.py`.
    - The header's "Model Running" tab and most form controls are
      placeholders (no backend wired yet). They're inert — clicking them
      does nothing — rather than navigating or submitting data; see the
      `!Hidden` comments at each one for what they used to trigger.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpacerItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from client.workers import StartRunWorker, UploadWorker
from styles.theme import (
    LABEL_CAPS_STYLE,
    PAGE_MARGIN,
    PAGE_TITLE_STYLE,
    SECTION_SPACING,
    START_ACTION_BUTTON_STYLE,
    TEXT_MUTED,
)
from widgets.cards import SurfaceCard, UploadCard
from widgets.icons import icon_text
from widgets.navigation import (
    build_header_bar,
    build_sidebar,
    make_icon_button,
    make_primary_cta_button,
    make_sidebar_nav_button,
    make_top_tab,
)


class DatasetInitializationPage(QWidget):
    """Lets the user stage a dataset upload and configure an analysis run."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_initialize_upload: Callable[[dict], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_initialize_upload: Invoked with the backend's upload+run-start
                response (`run_id`, `device`, `expected_duration_seconds`,
                ...) once the CSV is uploaded and a model run has started;
                should navigate to the model execution log page. If `None`,
                the button shows a placeholder message instead.
        """
        super().__init__(parent)
        self.setObjectName("DatasetInitializationPage")
        self._on_initialize_upload = on_initialize_upload
        self._selected_file_path: str | None = None
        self._upload_heading: QLabel | None = None
        self._initialize_button: QPushButton | None = None
        self._upload_worker: UploadWorker | None = None
        self._start_run_worker: StartRunWorker | None = None
        self._cell_line_combo: QComboBox | None = None
        self._selected_cell_line: str | None = None
        # "upload" until a CSV has been accepted and its cell-line list
        # fetched; "start_run" once the user just needs to pick a target
        # sample and kick off the run.
        self._stage = "upload"
        self._build_ui()

    def _build_ui(self) -> None:
        """Lay out the sidebar, header, and scrollable two-column form body."""
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        root.addWidget(self._build_sidebar(), 0)

        content_shell = QFrame()
        content_shell.setObjectName("RootShell")
        content_shell_layout = QVBoxLayout(content_shell)
        content_shell_layout.setContentsMargins(0, 0, 0, 0)
        content_shell_layout.setSpacing(0)

        content_shell_layout.addWidget(self._build_header(), 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN)
        content_layout.setSpacing(SECTION_SPACING)

        title = QLabel("Dataset Initialization")
        title.setStyleSheet(PAGE_TITLE_STYLE)
        content_layout.addWidget(title, 0, Qt.AlignmentFlag.AlignLeft)

        content_layout.addLayout(self._build_form_stack())

        scroll.setWidget(content)
        content_shell_layout.addWidget(scroll, 1)
        root.addWidget(content_shell, 1)

    def _build_form_stack(self) -> QVBoxLayout:
        """Build the single full-width column of functional cards.

        A big upload dropzone on top, then the "Target Sample" card and the
        "Initialize Upload" button side by side beneath it — each card fills
        the page's width rather than leaving a sparse second column, since
        there are only two functional steps on this page.
        """
        stack = QVBoxLayout()
        stack.setSpacing(SECTION_SPACING)
        stack.addWidget(self._build_upload_card())

        bottom_row = QHBoxLayout()
        bottom_row.setSpacing(SECTION_SPACING)
        bottom_row.addWidget(self._build_target_sample_card(), 2)
        bottom_row.addWidget(self._build_initialize_button(), 1)
        stack.addLayout(bottom_row)

        stack.addItem(QSpacerItem(0, 0, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))
        return stack

    def _build_sidebar(self) -> QFrame:
        """Build the shared sidebar.

        Unlike the other pages, this sidebar has no "Model Visualization" /
        "Model Logs" nav section (no model run exists yet).
        """
        # !Hidden: was wired to a "frontend placeholder only" message box.
        cta = make_primary_cta_button()
        footer_widgets = [make_sidebar_nav_button("Support", "help_outline")]
        return build_sidebar(footer_widgets=footer_widgets, cta_widget=cta)

    def _build_header(self) -> QFrame:
        """Build the shared header bar with "Upload" as the active tab.

        The "Results" tab is a `QToolButton` with a dropdown menu (instead
        of a plain nav button) standing in for a future results-history
        picker; `build_header_bar` accepts any widget type in its tab list.
        """
        upload_tab = make_top_tab("Upload", active=True)
        model_tab = make_top_tab("Model Running")
        results_button = self._build_results_menu_button()

        # !Hidden: both were wired to "frontend placeholder" message boxes.
        notifications = make_icon_button("notifications", "Notifications")
        account = make_icon_button("account_circle", "Account")
        return build_header_bar([upload_tab, model_tab, results_button], [notifications, account])

    def _build_results_menu_button(self) -> QToolButton:
        """Build the header's "Results" dropdown, listing placeholder result sets."""
        results_button = QToolButton()
        results_button.setObjectName("ResultsButton")
        results_button.setText("Results  " + icon_text("keyboard_arrow_down"))
        results_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)

        menu = QMenu(results_button)
        for label in ("Transcriptome—Proteome", "Transcriptome—Epigenome", "Proteome—"):
            action = QAction(label, results_button)
            # !Hidden: was wired to a "frontend placeholder" message box.
            menu.addAction(action)
        results_button.setMenu(menu)
        return results_button

    def _build_upload_card(self) -> UploadCard:
        """Build the dataset upload dropzone card."""
        card = UploadCard()
        layout = QVBoxLayout(card)
        layout.setContentsMargins(32, 32, 32, 32)
        layout.setSpacing(14)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        card.setMinimumHeight(360)

        icon = QLabel(icon_text("cloud_upload"))
        icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        icon.setStyleSheet("font-size: 56px; color: #747878; background: transparent; border: none;")

        heading = QLabel("Drag and drop files here")
        heading.setStyleSheet("font-size: 24px; font-weight: 600; color: #1a1c1c; background: transparent; border: none;")
        heading.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._upload_heading = heading

        browse = QPushButton("Browse Files")
        browse.setObjectName("SecondaryActionButton")
        browse.clicked.connect(self._on_browse_files)
        browse.setCursor(Qt.CursorShape.PointingHandCursor)

        card.fileDropped.connect(self._handle_file_selected)

        layout.addWidget(icon, 0, Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(heading, 0, Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(browse, 0, Qt.AlignmentFlag.AlignCenter)
        return card

    def _build_target_sample_card(self) -> SurfaceCard:
        """Build the "Target Sample" card: picks the cell line to run predictions for.

        Disabled until a dataset has been uploaded, since the choices come
        from the intersection of the uploaded CSV's cell lines and the
        model's fixed omics reference files (see `backend/routers/dataset.py`).
        """
        card = SurfaceCard()
        card.setMinimumHeight(140)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(12)

        title = QLabel("Target Sample")
        title.setStyleSheet(LABEL_CAPS_STYLE)
        subtitle = QLabel("The cell line predictions will be generated for, chosen from your uploaded dataset.")
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet(f"font-size: 13px; color: {TEXT_MUTED}; background: transparent; border: none;")
        layout.addWidget(title)
        layout.addWidget(subtitle)

        combo = QComboBox()
        combo.addItem("Upload a dataset first")
        combo.setEnabled(False)
        combo.setMinimumHeight(36)
        self._cell_line_combo = combo
        layout.addWidget(combo)
        layout.addStretch(1)
        return card

    def _build_initialize_button(self) -> QWidget:
        """Build the "Initialize Upload" submit button.

        Navigates via `on_initialize_upload` if supplied, otherwise shows a
        placeholder message.
        """
        container = QWidget()
        container.setMinimumHeight(140)
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        button = QPushButton("Initialize Upload")
        button.setObjectName("PrimaryActionButton")
        button.setCursor(Qt.CursorShape.PointingHandCursor)
        # Matches the Target Sample card's height so the two sit flush
        # side by side instead of this one floating short in the middle.
        button.setMinimumHeight(140)
        button.clicked.connect(self._on_initialize_clicked)
        self._initialize_button = button
        layout.addWidget(button)
        return container

    def _show_message(self, title: str, text: str) -> None:
        """Show a placeholder informational dialog for not-yet-wired controls."""
        QMessageBox.information(self, title, text)

    def _on_browse_files(self) -> None:
        """Open a file picker restricted to CSV files and stage the choice."""
        path, _ = QFileDialog.getOpenFileName(self, "Select Dataset CSV", "", "CSV Files (*.csv)")
        if path:
            self._handle_file_selected(path)

    def _handle_file_selected(self, path: str) -> None:
        """Stage `path` as the file to upload and reflect it in the dropzone."""
        self._selected_file_path = path
        if self._upload_heading is not None:
            self._upload_heading.setText(Path(path).name)

    def _on_initialize_clicked(self) -> None:
        """Dispatch the primary button's click to the current stage.

        Stage "upload": uploads the staged CSV and, on success, populates
        the Target Sample dropdown from the response's `cell_line_ids` and
        advances to stage "start_run". Stage "start_run": starts a model run
        for the selected cell line and navigates via `on_initialize_upload`.
        """
        if self._stage == "upload":
            self._start_upload()
        else:
            self._start_run()

    def _start_upload(self) -> None:
        if self._selected_file_path is None:
            self._show_message("Initialize Upload", "Select a dataset CSV first.")
            return

        if self._initialize_button is not None:
            self._initialize_button.setEnabled(False)
            self._initialize_button.setText("Uploading...")

        self._upload_worker = UploadWorker(self._selected_file_path, parent=self)
        self._upload_worker.succeeded.connect(self._on_upload_succeeded)
        self._upload_worker.failed.connect(self._on_upload_failed)
        self._upload_worker.start()

    def _on_upload_succeeded(self, result: dict) -> None:
        cell_line_ids = result.get("cell_line_ids") or []
        if not cell_line_ids:
            self._reset_to_upload_stage()
            QMessageBox.critical(
                self,
                "Upload Failed",
                "None of this dataset's cell lines have matching omics reference data, "
                "so no target sample can be selected.",
            )
            return

        if self._cell_line_combo is not None:
            self._cell_line_combo.clear()
            self._cell_line_combo.addItems(cell_line_ids)
            self._cell_line_combo.setEnabled(True)

        self._stage = "start_run"
        if self._initialize_button is not None:
            self._initialize_button.setEnabled(True)
            self._initialize_button.setText("Start Run")
            self._initialize_button.setStyleSheet(START_ACTION_BUTTON_STYLE)

    def _on_upload_failed(self, message: str) -> None:
        self._reset_to_upload_stage()
        QMessageBox.critical(self, "Upload Failed", message)

    def _start_run(self) -> None:
        if self._cell_line_combo is None or not self._cell_line_combo.isEnabled():
            self._show_message("Start Run", "Upload a dataset first.")
            return
        self._selected_cell_line = self._cell_line_combo.currentText()

        if self._initialize_button is not None:
            self._initialize_button.setEnabled(False)
            self._initialize_button.setText("Starting Run...")

        self._start_run_worker = StartRunWorker(self._selected_cell_line, parent=self)
        self._start_run_worker.succeeded.connect(self._on_start_run_succeeded)
        self._start_run_worker.failed.connect(self._on_start_run_failed)
        self._start_run_worker.start()

    def _on_start_run_succeeded(self, result: dict) -> None:
        self._reset_to_upload_stage()
        if self._on_initialize_upload is not None:
            self._on_initialize_upload({**result, "target_cell_line": self._selected_cell_line})
        else:
            self._show_message("Start Run", "Run started (no navigation callback configured).")

    def _on_start_run_failed(self, message: str) -> None:
        if self._initialize_button is not None:
            self._initialize_button.setEnabled(True)
            self._initialize_button.setText("Start Run")
        QMessageBox.critical(self, "Run Failed to Start", message)

    def _reset_to_upload_stage(self) -> None:
        self._stage = "upload"
        if self._cell_line_combo is not None:
            self._cell_line_combo.clear()
            self._cell_line_combo.addItem("Upload a dataset first")
            self._cell_line_combo.setEnabled(False)
        if self._initialize_button is not None:
            self._initialize_button.setEnabled(True)
            self._initialize_button.setText("Initialize Upload")
            self._initialize_button.setStyleSheet("")
