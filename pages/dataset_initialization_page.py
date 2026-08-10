"""Dataset Initialization page.

The first step in the workflow: lets a user stage a dataset upload, verify a
target drug, select required cancer/omics types, and configure an analysis
profile before kicking off a model run.

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) with an extra
      page-specific "settings shortcuts" group between the nav section and
      the footer links.
    - Shared header bar (`widgets.navigation.build_header_bar`) with
      "Upload" marked as the active tab and a "Results" dropdown menu
      (`QToolButton`) standing in for a future results-history picker.
    - An upload dropzone card, target drug verification card, cancer-type
      card, omics-type selector, and analysis profile card, arranged in a
      two-column grid.

Interactions with other pages:
    - `on_initialize_upload` navigates to `ModelExecutionLogPage`, wired by
      `main.py`.
    - The header's "Model Running" tab and most form controls are
      placeholders (no backend wired yet) — clicking them shows an
      informational message box via `_show_message` rather than navigating
      or submitting data.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import partial

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
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

from styles.theme import LABEL_CAPS_STYLE, PAGE_MARGIN, PAGE_TITLE_STYLE, SECTION_SPACING
from widgets.cards import OptionTile, SurfaceCard, UploadCard
from widgets.icons import icon_text
from widgets.navigation import (
    build_header_bar,
    build_sidebar,
    make_icon_button,
    make_primary_cta_button,
    make_sidebar_nav_button,
    make_top_tab,
)

# Placeholder shortcut labels shown in the sidebar's settings group.
_SETTINGS_SHORTCUTS = ("Setting 1", "Setting 2", "Setting 3")


class DatasetInitializationPage(QWidget):
    """Lets the user stage a dataset upload and configure an analysis run."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_initialize_upload: Callable[[], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_initialize_upload: Invoked when "Initialize Upload" is
                clicked; should navigate to the model execution log page.
                If `None`, the button shows a placeholder message instead.
        """
        super().__init__(parent)
        self.setObjectName("DatasetInitializationPage")
        self._omics_tiles: list[OptionTile] = []
        self._on_initialize_upload = on_initialize_upload
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

        content_layout.addLayout(self._build_form_grid())

        scroll.setWidget(content)
        content_shell_layout.addWidget(scroll, 1)
        root.addWidget(content_shell, 1)

    def _build_form_grid(self) -> QGridLayout:
        """Build the two-column grid of configuration cards.

        Left column: upload dropzone, target drug verification, cancer
        types. Right column: omics type selector, analysis profile, the
        "Initialize Upload" submit button.
        """
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(SECTION_SPACING)
        grid.setVerticalSpacing(SECTION_SPACING)
        grid.setColumnStretch(0, 8)
        grid.setColumnStretch(1, 4)

        left_stack = QVBoxLayout()
        left_stack.setSpacing(SECTION_SPACING)
        left_stack.addWidget(self._build_upload_card())
        left_stack.addWidget(self._build_target_verification_card())
        left_stack.addWidget(self._build_cancer_types_card())
        left_stack.addItem(QSpacerItem(0, 0, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        right_stack = QVBoxLayout()
        right_stack.setSpacing(SECTION_SPACING)
        right_stack.addWidget(self._build_omics_card())
        right_stack.addWidget(self._build_analysis_profile_card())
        right_stack.addWidget(self._build_initialize_button())
        right_stack.addItem(QSpacerItem(0, 0, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

        left_host = QWidget()
        left_host.setLayout(left_stack)
        right_host = QWidget()
        right_host.setLayout(right_stack)

        grid.addWidget(left_host, 0, 0)
        grid.addWidget(right_host, 0, 1)
        return grid

    def _build_sidebar(self) -> QFrame:
        """Build the shared sidebar with a page-specific settings shortcuts group.

        Unlike the other pages, this sidebar has no "Model Visualization" /
        "Model Logs" nav section (no model run exists yet); it instead shows
        a small settings-shortcut group above the standard footer links.
        """
        cta = make_primary_cta_button(
            callback=partial(self._show_message, "New Analysis", "This page uses a frontend placeholder only.")
        )
        footer_widgets = [
            make_sidebar_nav_button("History", "history"),
            make_sidebar_nav_button("Settings", "settings"),
            make_sidebar_nav_button("Support", "help_outline"),
        ]
        return build_sidebar(
            footer_widgets=footer_widgets,
            cta_widget=cta,
            extra_widget=self._build_settings_shortcuts(),
        )

    def _build_settings_shortcuts(self) -> QWidget:
        """Build the sidebar's page-specific settings shortcut links and separator."""
        group = QWidget()
        group_layout = QVBoxLayout(group)
        group_layout.setContentsMargins(0, 0, 0, 0)
        group_layout.setSpacing(6)

        for label in _SETTINGS_SHORTCUTS:
            group_layout.addWidget(make_sidebar_nav_button(label, "settings"))

        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setObjectName("SidebarSeparator")
        group_layout.addWidget(separator)
        return group

    def _build_header(self) -> QFrame:
        """Build the shared header bar with "Upload" as the active tab.

        The "Results" tab is a `QToolButton` with a dropdown menu (instead
        of a plain nav button) standing in for a future results-history
        picker; `build_header_bar` accepts any widget type in its tab list.
        """
        upload_tab = make_top_tab("Upload", active=True)
        model_tab = make_top_tab("Model Running")
        results_button = self._build_results_menu_button()

        notifications = make_icon_button(
            "notifications", "Notifications", callback=partial(self._show_message, "Notifications", "Notifications is a frontend placeholder.")
        )
        account = make_icon_button(
            "account_circle", "Account", callback=partial(self._show_message, "Account", "Account is a frontend placeholder.")
        )
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
            action.triggered.connect(partial(self._show_message, "Results", f"{label} is a frontend placeholder."))
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
        card.setMinimumHeight(300)

        icon = QLabel(icon_text("cloud_upload"))
        icon.setAlignment(Qt.AlignmentFlag.AlignCenter)
        icon.setStyleSheet("font-size: 56px; color: #747878; background: transparent; border: none;")

        heading = QLabel("Drag and drop files here")
        heading.setStyleSheet("font-size: 24px; font-weight: 600; color: #1a1c1c; background: transparent; border: none;")
        heading.setAlignment(Qt.AlignmentFlag.AlignCenter)

        browse = QPushButton("Browse Files")
        browse.setObjectName("SecondaryActionButton")
        browse.clicked.connect(partial(self._show_message, "Browse Files", "File browsing is a frontend placeholder."))
        browse.setCursor(Qt.CursorShape.PointingHandCursor)

        layout.addWidget(icon, 0, Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(heading, 0, Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(browse, 0, Qt.AlignmentFlag.AlignCenter)
        return card

    def _build_target_verification_card(self) -> SurfaceCard:
        """Build the target drug verification card (search box + select-all)."""
        card = SurfaceCard()
        layout = QVBoxLayout(card)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(16)

        header = QHBoxLayout()
        title = QLabel("Target Drug Verification")
        title.setStyleSheet(LABEL_CAPS_STYLE)
        select_all = QCheckBox("Select All")
        header.addWidget(title)
        header.addStretch(1)
        header.addWidget(select_all)
        layout.addLayout(header)

        form = QVBoxLayout()
        form.setSpacing(12)

        drug_label = QLabel("Drug Name")
        drug_label.setStyleSheet(LABEL_CAPS_STYLE)
        drug_input_row = QHBoxLayout()
        drug_input_row.setSpacing(8)

        drug_input = self._make_line_edit("Type drug name to verify...")
        search_button = make_icon_button(
            "search", "Search", callback=partial(self._show_message, "Drug Search", "Drug verification is a frontend placeholder.")
        )

        drug_input_row.addWidget(drug_input, 1)
        drug_input_row.addWidget(search_button, 0, Qt.AlignmentFlag.AlignVCenter)

        form.addWidget(drug_label)
        form.addLayout(drug_input_row)
        layout.addLayout(form)
        return card

    def _build_cancer_types_card(self) -> SurfaceCard:
        """Build the required cancer types input card."""
        card = SurfaceCard()
        layout = QVBoxLayout(card)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(12)

        title = QLabel("Required Cancer Types")
        title.setStyleSheet(LABEL_CAPS_STYLE)
        cancer_input = self._make_line_edit("Type cancer types...")

        layout.addWidget(title)
        layout.addWidget(cancer_input)
        return card

    def _build_omics_card(self) -> SurfaceCard:
        """Build the required omics types selector card.

        Populates `self._omics_tiles` so `_select_omics_tile` can enforce
        single-selection across the tiles.
        """
        card = SurfaceCard()
        layout = QVBoxLayout(card)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(12)

        title = QLabel("Required Omics Types")
        title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(title)

        for index, label in enumerate(("Omic Matrix 1", "Omic Matrix 2", "Omic Matrix 3")):
            tile = OptionTile(label, checked=index == 0)
            tile.clicked.connect(partial(self._select_omics_tile, tile))
            self._omics_tiles.append(tile)
            layout.addWidget(tile)
        return card

    def _build_analysis_profile_card(self) -> SurfaceCard:
        """Build the analysis profile card (cohort ID + reference genome)."""
        card = SurfaceCard()
        layout = QVBoxLayout(card)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(12)

        title = QLabel("Analysis Profile")
        title.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(title)

        cohort_label = QLabel("Cohort ID")
        cohort_label.setStyleSheet(LABEL_CAPS_STYLE)
        cohort_input = QLineEdit()
        cohort_input.setPlaceholderText("e.g. COH-2023-A")
        cohort_input.setStyleSheet("font-size: 13px; font-family: 'Consolas', 'Cascadia Mono', monospace;")

        genome_label = QLabel("Reference Genome")
        genome_label.setStyleSheet(LABEL_CAPS_STYLE)
        genome_combo = QComboBox()
        genome_combo.addItems(["GRCh38 (hg38)", "GRCh37 (hg19)"])

        layout.addWidget(cohort_label)
        layout.addWidget(cohort_input)
        layout.addWidget(genome_label)
        layout.addWidget(genome_combo)
        return card

    def _build_initialize_button(self) -> QWidget:
        """Build the "Initialize Upload" submit button.

        Navigates via `on_initialize_upload` if supplied, otherwise shows a
        placeholder message.
        """
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        button = QPushButton("Initialize Upload")
        button.setObjectName("PrimaryActionButton")
        button.setCursor(Qt.CursorShape.PointingHandCursor)
        if self._on_initialize_upload is None:
            button.clicked.connect(
                partial(self._show_message, "Initialize Upload", "Initialization is a frontend placeholder.")
            )
        else:
            button.clicked.connect(self._on_initialize_upload)
        layout.addWidget(button)
        return container

    @staticmethod
    def _make_line_edit(placeholder: str) -> QLineEdit:
        """Build a `QLineEdit` with the given placeholder text."""
        line_edit = QLineEdit()
        line_edit.setPlaceholderText(placeholder)
        return line_edit

    def _select_omics_tile(self, selected_tile: OptionTile) -> None:
        """Mark `selected_tile` as the sole selected omics tile.

        Args:
            selected_tile: The tile the user just clicked.
        """
        for tile in self._omics_tiles:
            tile.setSelected(tile is selected_tile)

    def _show_message(self, title: str, text: str) -> None:
        """Show a placeholder informational dialog for not-yet-wired controls."""
        QMessageBox.information(self, title, text)
