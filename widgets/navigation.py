"""Shared application "chrome": the sidebar and top header bar.

Every page in `pages/` shows the same brand-marked sidebar (logo, primary
call-to-action, navigation links, footer links) and the same top header bar
(workflow tabs plus notification/account icon buttons). Before this module
existed, each page file re-implemented that chrome with its own copy of the
layout and inline styles, which is how the pages drifted out of sync (e.g.
`final_results_page.py` grew a differently colored header).

Pages now build the *content* that is specific to them (which nav items are
shown, which tab is active, which callbacks fire) and hand it to the builder
functions here, which are responsible for the actual layout and styling.
That keeps the logo, spacing, and colors identical across the app by
construction rather than by convention.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget,
)

from styles.theme import (
    BUTTON_RADIUS,
    FONT_LABEL_CAPS,
    HEADER_HEIGHT,
    PAGE_MARGIN,
    PRIMARY,
    SIDEBAR_WIDTH,
    SURFACE,
    SURFACE_HIGH,
    TEXT,
    TEXT_MUTED,
)
from widgets.icons import icon_text


def build_brand_row() -> QHBoxLayout:
    """Build the "MSC16 / Biomedical Engine v2.4" logo block.

    Used at the top of every sidebar so the brand mark sits in the exact
    same position, with the same icon, spacing, and typography, on every
    page in the application.

    Returns:
        A `QHBoxLayout` containing the brand icon and title/subtitle column,
        ready to be added to a sidebar layout.
    """
    row = QHBoxLayout()
    row.setSpacing(12)

    icon = QLabel(icon_text("biotech"))
    icon.setStyleSheet(f"font-size: 32px; color: {TEXT}; background: transparent; border: none;")

    text_column = QVBoxLayout()
    text_column.setSpacing(2)
    title = QLabel("MSC16")
    title.setObjectName("BrandTitle")
    subtitle = QLabel("Biomedical Engine v2.4")
    subtitle.setStyleSheet(
        f"font-size: 12px; color: {TEXT_MUTED}; font-family: Consolas, monospace; background: transparent; border: none;"
    )
    text_column.addWidget(title)
    text_column.addWidget(subtitle)

    row.addWidget(icon, 0, Qt.AlignmentFlag.AlignTop)
    row.addLayout(text_column)
    row.addStretch(1)
    return row


def make_primary_cta_button(
    text: str = "New Analysis",
    icon_name: str = "add",
    callback: Callable[[], None] | None = None,
) -> QPushButton:
    """Build the compact, filled call-to-action button shown in the sidebar.

    Args:
        text: Button label.
        icon_name: Key into `widgets.icons.ICON_MAP` for the leading glyph.
        callback: Optional slot connected to `clicked`.

    Returns:
        A styled `QPushButton` ready to add to a sidebar layout.
    """
    button = QPushButton(f"{icon_text(icon_name)}  {text}")
    button.setCursor(Qt.CursorShape.PointingHandCursor)
    button.setStyleSheet(
        f"background: {PRIMARY}; color: {SURFACE}; border: none; border-radius: {BUTTON_RADIUS}px;"
        " padding: 8px 16px; font-size: 12px; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase;"
    )
    if callback is not None:
        button.clicked.connect(callback)
    return button


def make_sidebar_nav_button(
    text: str,
    icon_name: str,
    *,
    active: bool = False,
    callback: Callable[[], None] | None = None,
) -> QPushButton:
    """Build one sidebar navigation link (used for both nav and footer links).

    Args:
        text: Link label.
        icon_name: Key into `widgets.icons.ICON_MAP` for the leading glyph.
        active: Whether this link represents the currently displayed page.
            Active links get a filled background and primary text color.
        callback: Optional slot connected to `clicked`.

    Returns:
        A styled `QPushButton` ready to add to a sidebar layout.
    """
    background = SURFACE_HIGH if active else "transparent"
    color = PRIMARY if active else TEXT_MUTED

    button = QPushButton(f"{icon_text(icon_name)}  {text}")
    button.setCursor(Qt.CursorShape.PointingHandCursor)
    button.setStyleSheet(
        f"background: {background}; border: none; text-align: left; padding: 8px 12px; color: {color};"
        f" border-radius: {BUTTON_RADIUS}px; font-size: {FONT_LABEL_CAPS}px; font-weight: 700;"
        " letter-spacing: 0.05em; text-transform: uppercase;"
    )
    if callback is not None:
        button.clicked.connect(callback)
    return button


def build_sidebar(
    nav_widgets: Sequence[QWidget] = (),
    footer_widgets: Sequence[QWidget] = (),
    *,
    cta_widget: QWidget | None = None,
    extra_widget: QWidget | None = None,
) -> QFrame:
    """Assemble the full sidebar shell shared by every page.

    Layout order (top to bottom) is fixed for all pages: brand row, optional
    CTA button, optional primary navigation links, optional page-specific
    extra content, a flexible spacer, then footer links. Passing empty
    sequences/`None` for a section simply omits it, so pages that don't have
    a "Model Visualization / Model Logs" nav (e.g. the upload and results
    pages) still get an identically positioned logo and footer.

    Args:
        nav_widgets: Pre-built nav buttons (see `make_sidebar_nav_button`).
        footer_widgets: Pre-built footer links, rendered just above the
            bottom edge.
        cta_widget: Optional call-to-action button placed under the brand
            row (see `make_primary_cta_button`).
        extra_widget: Optional page-specific content placed after the nav
            section (e.g. the dataset page's settings shortcuts).

    Returns:
        A fully laid-out `QFrame` with object name "Sidebar", ready to be
        added to a page's root layout.
    """
    sidebar = QFrame()
    sidebar.setObjectName("Sidebar")
    sidebar.setFixedWidth(SIDEBAR_WIDTH)
    sidebar.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)

    layout = QVBoxLayout(sidebar)
    layout.setContentsMargins(24, 24, 24, 24)
    layout.setSpacing(16)

    layout.addLayout(build_brand_row())

    if cta_widget is not None:
        layout.addWidget(cta_widget)

    if nav_widgets:
        nav_layout = QVBoxLayout()
        nav_layout.setSpacing(6)
        for widget in nav_widgets:
            nav_layout.addWidget(widget)
        layout.addLayout(nav_layout)

    if extra_widget is not None:
        layout.addWidget(extra_widget)

    # Pushes the footer links to the bottom regardless of how much nav
    # content precedes them.
    layout.addItem(QSpacerItem(0, 0, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))

    if footer_widgets:
        footer_layout = QVBoxLayout()
        footer_layout.setSpacing(6)
        for widget in footer_widgets:
            footer_layout.addWidget(widget)
        layout.addLayout(footer_layout)

    return sidebar


def make_top_tab(
    text: str,
    *,
    active: bool = False,
    callback: Callable[[], None] | None = None,
) -> QPushButton:
    """Build one workflow tab ("Upload" / "Model Running" / "Results") for the header.

    Args:
        text: Tab label.
        active: Whether this tab represents the currently displayed page.
            The active tab is underlined and uses the primary text color.
        callback: Optional slot connected to `clicked`.

    Returns:
        A styled `QPushButton` ready to add to a header layout. Callers that
        need a dropdown menu instead of a plain click (see the dataset page's
        "Results" tab) can build a `QToolButton` directly and pass it to
        `build_header_bar` alongside these — the header layout doesn't care
        about the concrete widget type.
    """
    button = QPushButton(text)
    button.setCursor(Qt.CursorShape.PointingHandCursor)
    if active:
        button.setStyleSheet(
            f"background: transparent; border: none; color: {PRIMARY}; border-bottom: 2px solid {PRIMARY};"
            " padding: 0 8px 6px 8px; font-size: 12px; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase;"
        )
    else:
        button.setStyleSheet(
            f"background: transparent; border: none; color: {TEXT_MUTED}; padding: 0 8px 6px 8px;"
            " font-size: 12px; font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase;"
        )
    if callback is not None:
        button.clicked.connect(callback)
    return button


def make_icon_button(
    icon_name: str,
    tooltip: str = "",
    callback: Callable[[], None] | None = None,
) -> QPushButton:
    """Build a borderless icon-only button for header/toolbar actions.

    Reuses the global `QPushButton#IconButton` rule from `APP_STYLESHEET`
    (transparent by default, subtle hover background) so every icon button
    in the app — notifications, account, download, search — behaves and
    looks identical.

    Args:
        icon_name: Key into `widgets.icons.ICON_MAP`.
        tooltip: Optional hover tooltip text.
        callback: Optional slot connected to `clicked`.

    Returns:
        A styled `QPushButton`.
    """
    button = QPushButton(icon_text(icon_name))
    button.setObjectName("IconButton")
    button.setCursor(Qt.CursorShape.PointingHandCursor)
    if tooltip:
        button.setToolTip(tooltip)
    if callback is not None:
        button.clicked.connect(callback)
    return button


def build_header_bar(
    tab_widgets: Sequence[QWidget],
    trailing_widgets: Sequence[QWidget] | None = None,
) -> QFrame:
    """Assemble the top header bar shared by every page.

    Args:
        tab_widgets: Pre-built workflow tabs, in left-to-right order (see
            `make_top_tab`). Any widget type is accepted, so a page can mix
            in a `QToolButton` with a dropdown menu where needed.
        trailing_widgets: Pre-built widgets shown right-aligned, typically
            icon buttons. Defaults to the standard notifications + account
            pair via `make_icon_button` when omitted.

    Returns:
        A fully laid-out `QFrame` with object name "HeaderBar", ready to be
        added to a page's root layout.
    """
    header = QFrame()
    header.setObjectName("HeaderBar")
    header.setFixedHeight(HEADER_HEIGHT)
    header.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)

    layout = QHBoxLayout(header)
    layout.setContentsMargins(PAGE_MARGIN, 0, PAGE_MARGIN, 0)
    layout.setSpacing(16)

    tabs_layout = QHBoxLayout()
    tabs_layout.setSpacing(24)
    for tab in tab_widgets:
        tabs_layout.addWidget(tab)
    tabs_host = QWidget()
    tabs_host.setLayout(tabs_layout)
    layout.addWidget(tabs_host, 1)
    layout.addStretch(1)

    if trailing_widgets is None:
        trailing_widgets = (
            make_icon_button("notifications", "Notifications"),
            make_icon_button("account_circle", "Account"),
        )
    for widget in trailing_widgets:
        layout.addWidget(widget)

    return header
