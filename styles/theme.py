"""Application-wide visual theme.

This module is the single source of truth for the design language shared by
every page in `pages/`: colors, spacing, radii, and font sizes. Widgets in
`widgets/navigation.py` and `widgets/tables.py`, and the page modules
themselves, import these constants rather than hard-coding literals so that a
palette or spacing change only needs to happen in one place.

`apply_theme()` is called once from `main.py` to install the base Qt style
and the global stylesheet (`APP_STYLESHEET`) on the `QApplication` instance.
"""

from __future__ import annotations

from PySide6.QtGui import QFont
from PySide6.QtWidgets import QApplication

# --- Core palette -----------------------------------------------------------
# Neutral, near-monochrome palette shared by every page. Keep additions here
# rather than introducing new one-off hex literals in page modules.
WINDOW_BACKGROUND = "#f9f9f9"
SURFACE = "#ffffff"
SURFACE_LOW = "#f3f3f3"
SURFACE_CONTAINER = "#eeeeee"
SURFACE_HIGH = "#e8e8e8"
SURFACE_HIGHEST = "#e2e2e2"
BORDER = "#c4c7c7"
DIVIDER = "#e2e2e2"
TEXT = "#1a1c1c"
TEXT_MUTED = "#444748"
TEXT_SOFT = "#63646c"
TEXT_FAINT = "#a8abab"
PRIMARY = "#000000"
PRIMARY_SOFT = "#2f3131"
ACCENT = "#5f5e5e"
SUCCESS = "#1a7f37"
SUCCESS_SOFT = "#156429"

# --- Layout scale ------------------------------------------------------------
# Shared sizing so every page's chrome (sidebar/header) lines up pixel-for-pixel.
SIDEBAR_WIDTH = 260
HEADER_HEIGHT = 64
PAGE_MARGIN = 32
SECTION_SPACING = 24
CARD_PADDING = 24

# --- Radii --------------------------------------------------------------------
CARD_RADIUS = 12
BUTTON_RADIUS = 4
PILL_RADIUS = 999

# --- Typography scale ---------------------------------------------------------
FONT_PAGE_TITLE = 30
FONT_PAGE_SUBTITLE = 16
FONT_CARD_TITLE = 18
FONT_LABEL_CAPS = 12
FONT_BODY = 14
FONT_CODE = 13

# Reusable style fragments so pages compose consistent inline stylesheets
# instead of repeating the same declarations with slightly different values.
# Qt/Fusion quirk: a plain QLabel nested inside any ancestor widget that has
# its own inline stylesheet (most card headers/bodies do) can be painted
# with a phantom sunken border unless the label explicitly opts out with its
# own "background: transparent; border: none;". A global QLabel rule in
# APP_STYLESHEET does NOT reliably prevent this, because a styled ancestor's
# cascade takes precedence over an app-wide type selector — so every label
# style used across the app is built starting from this reset via
# `label_style()` rather than relying on the global rule alone.
LABEL_RESET_STYLE = "background: transparent; border: none;"


def label_style(extra: str) -> str:
    """Build a `QLabel` stylesheet that starts from the transparent/borderless reset.

    Args:
        extra: Additional QSS declarations (e.g. font/color rules) to apply
            after the reset.

    Returns:
        A complete stylesheet string for `QLabel.setStyleSheet()`.
    """
    return f"{LABEL_RESET_STYLE} {extra}"


PAGE_TITLE_STYLE = label_style(f"font-size: {FONT_PAGE_TITLE}px; font-weight: 600; color: {TEXT}; letter-spacing: -0.02em;")
PAGE_SUBTITLE_STYLE = label_style(f"font-size: {FONT_PAGE_SUBTITLE}px; color: {TEXT_MUTED};")
CARD_TITLE_STYLE = label_style(f"font-size: {FONT_CARD_TITLE}px; font-weight: 600; color: {TEXT};")
LABEL_CAPS_STYLE = label_style(f"font-size: {FONT_LABEL_CAPS}px; font-weight: 700; letter-spacing: 0.05em; color: {TEXT_MUTED};")
CARD_CONTAINER_STYLE = f"background: {SURFACE}; border: 1px solid {BORDER}; border-radius: {CARD_RADIUS}px;"

# Small, uppercase, fixed-height action buttons (e.g. "Pause Execution",
# "Halt Execution", "Download Data", "Export Clinical Report", "Finish").
# Sharing one height/padding/radius across pages keeps these visually
# identical regardless of which page renders them.
SECONDARY_BUTTON_STYLE = (
    f"height: 40px; padding: 0 16px; border: 1px solid {BORDER}; border-radius: {BUTTON_RADIUS}px;"
    f" background: {SURFACE}; color: {TEXT}; font-size: 12px; font-weight: 700; letter-spacing: 0.05em;"
    " text-transform: uppercase;"
)
PRIMARY_BUTTON_STYLE = (
    f"height: 40px; padding: 0 16px; border: none; border-radius: {BUTTON_RADIUS}px; background: {PRIMARY};"
    f" color: {SURFACE}; font-size: 12px; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase;"
)

# Green variant of the big `QPushButton#PrimaryActionButton` QSS rule (see
# APP_STYLESHEET below), set inline on a specific button instance to signal
# a "ready to go, action is affirmative" state — e.g. the dataset
# initialization page's submit button once it's about to start a run rather
# than upload a file.
START_ACTION_BUTTON_STYLE = f"""
    QPushButton {{
        background: {SUCCESS};
        color: white;
        border: none;
        border-radius: 10px;
        padding: 12px 16px;
        font-size: 18px;
        font-weight: 600;
    }}
    QPushButton:hover {{
        background: {SUCCESS_SOFT};
    }}
"""


def table_stylesheet(header_background: str = WINDOW_BACKGROUND) -> str:
    """Build the shared QSS used by every data table in the app.

    Centralizing this means the pipeline-run table and the predicted-drug
    results table (and any future table) render with identical header,
    border, and row styling regardless of which page hosts them.

    Args:
        header_background: Background color for the `QHeaderView` section.
            Defaults to the app window background; pass a different shared
            color (e.g. `SURFACE`) to match a specific card's header strip.

    Returns:
        A QSS string suitable for `QTableWidget.setStyleSheet()`.
    """
    return (
        f"QTableWidget {{ border: none; background: {SURFACE}; font-size: {FONT_CODE}px; }}"
        f"QHeaderView::section {{ background: {header_background}; color: {TEXT_MUTED}; padding: 12px 16px;"
        f" border: none; border-bottom: 1px solid {DIVIDER}; font-size: {FONT_LABEL_CAPS}px; font-weight: 700;"
        " letter-spacing: 0.05em; }"
        "QTableWidget::item { padding: 10px 16px; }"
        "QTableWidget::item:selected { background: transparent; }"
    )

APP_STYLESHEET = f"""
QWidget {{
    color: {TEXT};
    background: transparent;
    font-family: "Segoe UI", "Inter", sans-serif;
    font-size: 14px;
}}

/*
 * Qt/Fusion quirk: once any ancestor widget has its own inline stylesheet
 * (e.g. a card's `background`/`border`), a plain child QLabel with no
 * stylesheet of its own can be painted with a phantom sunken frame instead
 * of rendering as transparent text. Setting a global, low-priority default
 * here avoids having to repeat "background: transparent; border: none;" on
 * every single label — any label that DOES set its own background/border
 * (badges, pills, chips) simply overrides this via its own stylesheet.
 */
QLabel {{
    background: transparent;
    border: none;
}}

QMainWindow, QWidget#DatasetInitializationPage, QWidget#RootShell {{
    background: {WINDOW_BACKGROUND};
}}

QFrame#Sidebar {{
    background: {SURFACE};
    border-right: 1px solid {BORDER};
}}

QFrame#HeaderBar {{
    background: {SURFACE};
    border-bottom: 1px solid {BORDER};
}}

QFrame#SurfaceCard, QFrame#UploadCard, QFrame#OptionTile {{
    background: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 12px;
}}

QFrame#UploadCard {{
    border: 2px dashed {BORDER};
}}

QFrame#UploadCard:hover {{
    background: {SURFACE_LOW};
    border-color: {PRIMARY};
}}

QFrame#OptionTile[selected="true"] {{
    border: 2px solid {PRIMARY};
    background: rgba(0, 0, 0, 0.03);
}}

QFrame#OptionIndicator {{
    border: 2px solid {BORDER};
    border-radius: 10px;
    background: transparent;
}}

QFrame#OptionIndicator[selected="true"] {{
    border-color: {PRIMARY};
    background: {PRIMARY};
}}

QLabel#BrandTitle {{
    color: {PRIMARY};
    font-size: 24px;
    font-weight: 700;
    letter-spacing: -0.01em;
}}

QPushButton#PrimaryActionButton {{
    background: {PRIMARY};
    color: white;
    border: none;
    border-radius: 10px;
    padding: 12px 16px;
    font-size: 18px;
    font-weight: 600;
}}

QPushButton#PrimaryActionButton:hover {{
    background: {PRIMARY_SOFT};
}}

QPushButton#SecondaryActionButton {{
    background: {SURFACE};
    color: {PRIMARY};
    border: 1px solid {BORDER};
    border-radius: 8px;
    padding: 8px 16px;
    font-size: 18px;
    font-weight: 600;
}}

QPushButton#SecondaryActionButton:hover {{
    border-color: {PRIMARY};
    background: {SURFACE_LOW};
}}

QPushButton#IconButton {{
    background: transparent;
    color: {TEXT_SOFT};
    border: none;
    padding: 6px;
    border-radius: 8px;
    font-size: 20px;
}}

QPushButton#IconButton:hover {{
    color: {PRIMARY};
    background: {SURFACE_LOW};
}}

QToolButton#ResultsButton {{
    background: transparent;
    color: {TEXT_MUTED};
    border: none;
    padding: 0 8px 6px 8px;
    font-size: 12px;
    font-weight: 700;
    letter-spacing: 0.05em;
    text-transform: uppercase;
}}

QToolButton#ResultsButton::menu-indicator {{
    image: none;
}}

QLineEdit, QComboBox {{
    background: {SURFACE};
    color: {PRIMARY};
    border: 1px solid {BORDER};
    border-radius: 8px;
    padding: 8px 10px;
    selection-background-color: {PRIMARY};
    selection-color: white;
}}

QLineEdit:focus, QComboBox:focus {{
    border: 1px solid {PRIMARY};
}}

QComboBox::drop-down {{
    border: none;
    width: 24px;
}}

QComboBox QAbstractItemView {{
    background: {SURFACE};
    selection-background-color: {SURFACE_CONTAINER};
    outline: 0;
}}

QCheckBox {{
    color: {TEXT_SOFT};
    spacing: 8px;
}}

QCheckBox::indicator {{
    width: 14px;
    height: 14px;
    border-radius: 3px;
    border: 1px solid {BORDER};
    background: {SURFACE};
}}

QCheckBox::indicator:checked {{
    background: {PRIMARY};
    border-color: {PRIMARY};
}}

QRadioButton {{
    color: {TEXT};
    spacing: 8px;
}}

QRadioButton::indicator {{
    width: 18px;
    height: 18px;
    border-radius: 9px;
    border: 2px solid {BORDER};
    background: {SURFACE};
}}

QRadioButton::indicator:checked {{
    border-color: {PRIMARY};
    background: {PRIMARY};
}}

QMenu {{
    background: {SURFACE};
    border: 1px solid {BORDER};
    padding: 4px;
}}

QMenu::item {{
    padding: 8px 14px;
    border-radius: 6px;
}}

QMenu::item:selected {{
    background: {SURFACE_CONTAINER};
    color: {PRIMARY};
}}

QScrollBar:vertical {{
    width: 10px;
    background: transparent;
    margin: 0;
}}

QScrollBar::handle:vertical {{
    background: rgba(0, 0, 0, 0.18);
    min-height: 24px;
    border-radius: 5px;
}}

QScrollBar::handle:vertical:hover {{
    background: rgba(0, 0, 0, 0.28);
}}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}

QMessageBox {{
    background: {SURFACE};
}}
"""


def apply_theme(app: QApplication) -> None:
    """Install the app-wide Qt style, base font, and global stylesheet.

    Call this once, immediately after constructing the `QApplication`, so
    every widget created afterwards (sidebars, headers, cards, tables across
    all pages) picks up the same Fusion style and `APP_STYLESHEET` rules.

    Args:
        app: The application instance to theme.
    """
    app.setStyle("Fusion")
    app.setFont(QFont("Segoe UI", 10))
    app.setStyleSheet(APP_STYLESHEET)
