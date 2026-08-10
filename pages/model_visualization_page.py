"""Model Visualization page.

Shown while a model run is in progress, alongside the model execution log
page. It renders an animated "constellation" canvas standing in for the
model's growing computation graph, plus summary panels describing the run's
input source and architecture.

Main UI components:
    - Shared sidebar (`widgets.navigation.build_sidebar`) with "Model
      Visualization" marked as the active nav item.
    - Shared header bar (`widgets.navigation.build_header_bar`) with "Model
      Running" marked as the active workflow tab (visualization and logs
      are both views onto the same running-model stage).
    - `ConstellationCanvas`: a custom-painted, self-animating widget that
      spawns and connects nodes over time to suggest a growing network.
    - Two summary panels ("Input Source", "Architecture") built from
      key/value rows via `_kv`.

Interactions with other pages:
    - `on_upload_clicked` navigates back to `DatasetInitializationPage`.
    - `on_model_running_clicked` and the sidebar's "Model Logs" link both
      navigate to `ModelExecutionLogPage`.
    - `on_finish_clicked` navigates to `FinalResultsPage`, triggered either
      by the "Finish" button or automatically when `ConstellationCanvas`
      finishes spawning nodes (`completed` signal).
    All callbacks are supplied and wired by `main.py`.
"""

from __future__ import annotations

import math
import random
from collections.abc import Callable

from PySide6.QtCore import QElapsedTimer, QPointF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from styles.theme import (
    BORDER,
    CARD_CONTAINER_STYLE,
    LABEL_CAPS_STYLE,
    PAGE_MARGIN,
    PAGE_SUBTITLE_STYLE,
    PAGE_TITLE_STYLE,
    PRIMARY,
    PRIMARY_BUTTON_STYLE,
    SECONDARY_BUTTON_STYLE,
    SECTION_SPACING,
    SURFACE_CONTAINER,
    TEXT,
    WINDOW_BACKGROUND,
    label_style,
)
from widgets.icons import icon_text
from widgets.navigation import (
    build_header_bar,
    build_sidebar,
    make_primary_cta_button,
    make_sidebar_nav_button,
    make_top_tab,
)


class ConstellationCanvas(QWidget):
    """Animated placeholder visualization standing in for a live model/graph feed.

    Periodically spawns nodes at random positions around a central hub,
    connects each to its nearest under-connected neighbor, and fades/pulses
    them in. Once `_max_nodes` is reached, spawning stops and `completed` is
    emitted so the host page can advance the workflow automatically.

    All the numeric tuning parameters (spawn interval, node radii, timings)
    are placeholder visual design choices, not derived from real model or
    graph data.
    """

    completed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        # * This animation is a placeholder visual; replace with real graph/model-state feed later.
        self.setMinimumSize(700, 460)

        self._svg_width = 800
        self._svg_height = 600
        self._hub_x = 400
        self._hub_y = 300
        self._safe_margin = 40
        self._min_distance = 42
        self._max_nodes = 45
        self._zoom_threshold = 38
        self._small_node_radius = 8.0
        self._large_node_radius = 13.5
        self._special_node_probability = 0.12
        self._generation_interval_ms = 1550
        self._fade_in_ms = 900
        self._edge_draw_ms = 1800
        self._spawn_pulse_ms = 1100
        self._breath_period_ms = 2600

        self._clock = QElapsedTimer()
        self._clock.start()

        self._nodes: list[dict[str, float | int | bool]] = [
            {
                "x": self._hub_x,
                "y": self._hub_y,
                "connections": 0,
                "radius": self._small_node_radius,
                "special": False,
                "born_ms": 0,
            }
        ]
        self._edges: list[dict[str, int]] = []
        self._scale = 1.0

        # Drives node spawning; stops itself once `_max_nodes` is reached.
        self._timer = QTimer(self)
        self._timer.setInterval(self._generation_interval_ms)
        self._timer.timeout.connect(self._add_node)
        self._timer.start()

        # Drives repaint-only animation (fades, pulses, edge draw-in) at ~30fps.
        self._frame_timer = QTimer(self)
        self._frame_timer.setInterval(33)
        self._frame_timer.timeout.connect(self.update)
        self._frame_timer.start()

        self._add_node()

    @staticmethod
    def _clamp(value: float, min_value: float, max_value: float) -> float:
        """Clamp `value` into the inclusive range [min_value, max_value]."""
        return max(min_value, min(max_value, value))

    def _distance(self, x1: float, y1: float, x2: float, y2: float) -> float:
        """Return the Euclidean distance between two points."""
        return math.hypot(x2 - x1, y2 - y1)

    def _get_random_position(self) -> tuple[float, float] | None:
        """Sample a random position around the hub that isn't too close to an existing node.

        Returns:
            An (x, y) tuple in canvas coordinates, or `None` if no
            sufficiently spaced position was found after 120 attempts (the
            canvas is considered "full" for this spawn tick).
        """
        for _ in range(120):
            angle = random.random() * math.pi * 2.0
            radius = 50 + random.random() * (min(self._svg_width, self._svg_height) / 2 - self._safe_margin)

            x = self._hub_x + radius * math.cos(angle)
            y = self._hub_y + radius * math.sin(angle)

            x = max(self._safe_margin, min(self._svg_width - self._safe_margin, x))
            y = max(self._safe_margin, min(self._svg_height - self._safe_margin, y))

            if all(
                self._distance(x, y, float(node["x"]), float(node["y"])) >= self._min_distance
                for node in self._nodes
            ):
                return x, y
        return None

    def _find_parent_index(self, new_x: float, new_y: float) -> int:
        """Find the index of the nearest existing node with fewer than 2 connections.

        New nodes attach to this "parent" so the graph reads as a branching
        tree rather than a random scatter.
        """
        closest_index = 0
        min_distance = float("inf")
        for index, node in enumerate(self._nodes):
            if int(node["connections"]) < 2:
                distance = self._distance(new_x, new_y, float(node["x"]), float(node["y"]))
                if distance < min_distance:
                    min_distance = distance
                    closest_index = index
        return closest_index

    def _add_node(self) -> None:
        """Spawn one new node connected to its nearest under-connected parent.

        Stops the spawn timer and emits `completed` once `_max_nodes` is
        reached. Side effect: triggers a repaint via `update()`.
        """
        if len(self._nodes) >= self._max_nodes:
            self._timer.stop()
            self.completed.emit()
            return

        position = self._get_random_position()
        if position is None:
            return

        x, y = position
        parent_index = self._find_parent_index(x, y)

        is_special = random.random() < self._special_node_probability
        node_radius = self._large_node_radius if is_special else self._small_node_radius
        now_ms = int(self._clock.elapsed())

        self._nodes[parent_index]["connections"] = int(self._nodes[parent_index]["connections"]) + 1
        self._nodes.append(
            {
                "x": x,
                "y": y,
                "connections": 1,
                "radius": node_radius,
                "special": is_special,
                "born_ms": now_ms,
            }
        )
        self._edges.append({"parent": parent_index, "child": len(self._nodes) - 1, "born_ms": now_ms})

        # Zoom out gradually once the graph gets dense, to keep it in frame.
        if len(self._nodes) >= self._zoom_threshold:
            self._scale = max(0.6, 1 - ((len(self._nodes) - self._zoom_threshold) * 0.05))

        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        """Paint the dot-grid background, fading/growing edges, and pulsing nodes."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        painter.fillRect(self.rect(), QColor(WINDOW_BACKGROUND))

        dot_pen = QPen(QColor(212, 212, 216, 150))
        dot_pen.setWidth(1)
        painter.setPen(dot_pen)
        spacing = 24
        for y in range(0, self.height(), spacing):
            for x in range(0, self.width(), spacing):
                painter.drawPoint(x, y)

        # Fit the fixed-aspect-ratio SVG-space canvas into the widget, letterboxing as needed.
        canvas_w = self.width()
        canvas_h = self.height()
        target_ratio = self._svg_width / self._svg_height
        canvas_ratio = canvas_w / max(1, canvas_h)

        if canvas_ratio > target_ratio:
            draw_h = canvas_h
            draw_w = int(draw_h * target_ratio)
        else:
            draw_w = canvas_w
            draw_h = int(draw_w / target_ratio)

        offset_x = (canvas_w - draw_w) / 2
        offset_y = (canvas_h - draw_h) / 2

        painter.save()
        painter.translate(offset_x, offset_y)
        sx = draw_w / self._svg_width
        sy = draw_h / self._svg_height
        painter.scale(sx, sy)

        painter.translate(self._hub_x, self._hub_y)
        painter.scale(self._scale, self._scale)
        painter.translate(-self._hub_x, -self._hub_y)

        now_ms = int(self._clock.elapsed())

        edge_pen = QPen(QColor("#747878"))
        edge_pen.setWidthF(2.2)
        edge_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        for edge in self._edges:
            parent = self._nodes[edge["parent"]]
            child = self._nodes[edge["child"]]

            # Edges draw themselves in over `_edge_draw_ms` from parent to child.
            edge_age = max(0, now_ms - int(edge["born_ms"]))
            draw_ratio = self._clamp(edge_age / self._edge_draw_ms, 0.0, 1.0)
            edge_alpha = int(100 + (110 * draw_ratio))
            edge_pen.setColor(QColor(116, 120, 120, edge_alpha))
            painter.setPen(edge_pen)

            parent_x = float(parent["x"])
            parent_y = float(parent["y"])
            child_x = float(child["x"])
            child_y = float(child["y"])
            draw_x = parent_x + ((child_x - parent_x) * draw_ratio)
            draw_y = parent_y + ((child_y - parent_y) * draw_ratio)
            painter.drawLine(
                QPointF(parent_x, parent_y),
                QPointF(draw_x, draw_y),
            )

        painter.setPen(Qt.PenStyle.NoPen)
        for node in self._nodes:
            node_age = max(0, now_ms - int(node["born_ms"]))
            fade = self._clamp(node_age / self._fade_in_ms, 0.0, 1.0)

            color = QColor("#4ade80") if bool(node["special"]) else QColor(PRIMARY)
            color.setAlpha(int(255 * fade))
            painter.setBrush(color)

            radius = float(node["radius"])
            if bool(node["special"]):
                # Special nodes gently "breathe" (radius oscillates) to draw attention.
                breath = math.sin((node_age / self._breath_period_ms) * math.pi * 2.0)
                radius += 0.9 * breath

            x = float(node["x"])
            y = float(node["y"])

            # Pulse halo for newly spawned nodes.
            pulse_ratio = self._clamp(node_age / self._spawn_pulse_ms, 0.0, 1.0)
            if pulse_ratio < 1.0 and int(node["born_ms"]) > 0:
                halo_alpha = int((1.0 - pulse_ratio) * 140)
                halo_pen = QPen(QColor(0, 0, 0, halo_alpha))
                halo_pen.setWidthF(1.4)
                painter.setPen(halo_pen)
                painter.setBrush(Qt.BrushStyle.NoBrush)
                halo_radius = radius + (pulse_ratio * 8.0)
                painter.drawEllipse(QPointF(x, y), halo_radius, halo_radius)
                painter.setPen(Qt.PenStyle.NoPen)

            painter.drawEllipse(QPointF(x, y), radius, radius)

        painter.restore()


class ModelVisualizationPage(QWidget):
    """Live, animated visualization of an in-progress model run."""

    def __init__(
        self,
        parent: QWidget | None = None,
        on_upload_clicked: Callable[[], None] | None = None,
        on_model_running_clicked: Callable[[], None] | None = None,
        on_model_analytics_clicked: Callable[[], None] | None = None,
        on_finish_clicked: Callable[[], None] | None = None,
    ) -> None:
        """Build the page.

        Args:
            parent: Optional Qt parent widget.
            on_upload_clicked: Invoked when the header's "Upload" tab is
                clicked; should navigate back to the dataset upload page.
            on_model_running_clicked: Invoked when the header's "Model
                Running" tab is clicked; should navigate to the model
                execution log page.
            on_model_analytics_clicked: Invoked when the sidebar's "Model
                Logs" link is clicked; should navigate to the model
                execution log page.
            on_finish_clicked: Invoked when the run visualization completes
                (or the user clicks "Finish"); should navigate to the final
                results page.
        """
        super().__init__(parent)
        self.setObjectName("ModelVisualizationPage")
        # * Finish transition is placeholder behavior until real completion criteria is available.
        self._on_upload_clicked = on_upload_clicked
        self._on_model_running_clicked = on_model_running_clicked
        self._on_model_analytics_clicked = on_model_analytics_clicked
        self._on_finish_clicked = on_finish_clicked
        self._build_ui()

    def _build_ui(self) -> None:
        """Lay out the sidebar, header, canvas card, and summary panel grid."""
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        root.addWidget(self._build_sidebar(), 0)

        content_shell = QFrame()
        content_shell.setObjectName("RootShell")
        content_layout = QVBoxLayout(content_shell)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(0)

        content_layout.addWidget(self._build_header(), 0)

        main = QWidget()
        main_layout = QVBoxLayout(main)
        main_layout.setContentsMargins(PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN, PAGE_MARGIN)
        main_layout.setSpacing(SECTION_SPACING)

        main_layout.addLayout(self._build_page_header())
        main_layout.addWidget(self._build_visualization_card(), 1)

        bottom_grid = QGridLayout()
        bottom_grid.setContentsMargins(0, 0, 0, 0)
        bottom_grid.setHorizontalSpacing(SECTION_SPACING)
        bottom_grid.addWidget(self._build_info_panel("Input Source"), 0, 0)
        bottom_grid.addWidget(self._build_info_panel("Architecture"), 0, 1)
        main_layout.addLayout(bottom_grid)

        content_layout.addWidget(main, 1)
        root.addWidget(content_shell, 1)

    def _build_sidebar(self) -> QFrame:
        """Build the shared sidebar with "Model Visualization" as the active nav item."""
        nav_widgets = [
            make_sidebar_nav_button("Model Visualization", "hub", active=True),
            make_sidebar_nav_button("Model Logs", "terminal", callback=self._on_model_analytics_clicked),
        ]
        footer_widgets = [
            make_sidebar_nav_button("History", "history"),
            make_sidebar_nav_button("Settings", "settings"),
            make_sidebar_nav_button("Support", "help_outline"),
        ]
        return build_sidebar(nav_widgets, footer_widgets, cta_widget=make_primary_cta_button())

    def _build_header(self) -> QFrame:
        """Build the shared header bar with "Model Running" as the active tab."""
        tabs = [
            make_top_tab("Upload", callback=self._on_upload_clicked),
            make_top_tab("Model Running", active=True, callback=self._on_model_running_clicked),
            make_top_tab("Results"),
        ]
        return build_header_bar(tabs)

    def _build_page_header(self) -> QHBoxLayout:
        """Build the title/subtitle row plus the Finish/Halt action buttons."""
        title_row = QHBoxLayout()
        title_text = QVBoxLayout()
        title = QLabel("Execution Pipeline Active")
        title.setStyleSheet(PAGE_TITLE_STYLE)
        subtitle = QLabel("Processing cohort sequence data against baseline architecture.")
        subtitle.setStyleSheet(PAGE_SUBTITLE_STYLE)
        title_text.addWidget(title)
        title_text.addWidget(subtitle)
        title_row.addLayout(title_text)
        title_row.addStretch(1)

        finish_button = QPushButton("Finish")
        finish_button.setStyleSheet(SECONDARY_BUTTON_STYLE)
        if self._on_finish_clicked is not None:
            finish_button.clicked.connect(self._on_finish_clicked)

        halt_button = QPushButton(f"{icon_text('stop_circle')}  Halt Execution")
        halt_button.setStyleSheet(PRIMARY_BUTTON_STYLE)

        title_row.addWidget(finish_button)
        title_row.addWidget(halt_button)
        return title_row

    def _build_visualization_card(self) -> QFrame:
        """Build the card hosting the constellation canvas and its progress bar."""
        card = QFrame()
        card.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        canvas_host = QFrame()
        canvas_host.setStyleSheet(f"background: {WINDOW_BACKGROUND}; border-top-left-radius: 12px; border-top-right-radius: 12px;")
        canvas_layout = QVBoxLayout(canvas_host)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        canvas_layout.setSpacing(0)

        canvas = ConstellationCanvas()
        if self._on_finish_clicked is not None:
            canvas.completed.connect(self._on_finish_clicked)
        canvas_layout.addWidget(canvas)
        layout.addWidget(canvas_host, 1)

        # Static placeholder progress bar representing overall run completion.
        progress = QFrame()
        progress.setFixedHeight(4)
        progress.setStyleSheet(f"background: {SURFACE_CONTAINER}; border: none;")
        progress_layout = QHBoxLayout(progress)
        progress_layout.setContentsMargins(0, 0, 0, 0)
        progress_layout.setSpacing(0)
        fill = QFrame()
        fill.setFixedHeight(4)
        fill.setStyleSheet(f"background: {PRIMARY};")
        progress_layout.addWidget(fill, 68)
        progress_layout.addStretch(32)
        layout.addWidget(progress, 0)

        return card

    def _build_info_panel(self, panel_name: str) -> QFrame:
        """Build a small key/value summary card.

        Args:
            panel_name: Either `"Input Source"` or `"Architecture"`; selects
                both the icon and the set of rows shown.

        Returns:
            A styled `QFrame` card.
        """
        panel = QFrame()
        panel.setStyleSheet(CARD_CONTAINER_STYLE)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        icon = icon_text("folder_data") if panel_name == "Input Source" else icon_text("architecture")
        heading = QLabel(f"{icon}  {panel_name}")
        heading.setStyleSheet(LABEL_CAPS_STYLE)
        layout.addWidget(heading)

        if panel_name == "Input Source":
            layout.addWidget(self._kv("Cohort", "PT-99A", boxed=True))
            layout.addWidget(self._kv("Sample Size", "14,204"))
            layout.addWidget(self._kv("Sequence Type", "scRNA-seq", pill=True))
        else:
            layout.addWidget(self._kv("Base Model", "ResNet-Bio-V4"))
            layout.addWidget(self._kv("Epochs", "500"))
            layout.addWidget(self._kv("Learning Rate", "0.001"))
        return panel

    @staticmethod
    def _kv(key: str, value: str, boxed: bool = False, pill: bool = False) -> QWidget:
        """Build one key/value row for an info panel.

        Args:
            key: Row label, rendered in the shared label-caps style.
            value: Row value, rendered in monospace.
            boxed: Render the value with a filled rectangular background.
            pill: Render the value with a filled, fully-rounded background.
                Mutually exclusive with `boxed` in practice (callers pass at
                most one of the two).

        Returns:
            A `QWidget` wrapping the key/value column.
        """
        wrapper = QWidget()
        column = QVBoxLayout(wrapper)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(3)

        key_label = QLabel(key)
        key_label.setStyleSheet(LABEL_CAPS_STYLE)
        column.addWidget(key_label)

        value_label = QLabel(value)
        if boxed:
            value_label.setStyleSheet(
                label_style(f"font-family: Consolas, monospace; color: {TEXT}; background: {WINDOW_BACKGROUND};"
                            " border-radius: 4px; padding: 4px 8px;")
            )
        elif pill:
            value_label.setStyleSheet(
                label_style(f"font-family: Consolas, monospace; color: {TEXT}; background: {WINDOW_BACKGROUND};"
                            " border-radius: 12px; padding: 4px 10px;")
            )
        else:
            value_label.setStyleSheet(label_style(f"font-family: Consolas, monospace; color: {TEXT};"))
        column.addWidget(value_label)
        return wrapper
