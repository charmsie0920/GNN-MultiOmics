"""Application entry point for the MSC16 desktop frontend.

This is a PySide6 (Qt for Python) desktop application for a multi-omics
drug-response prediction tool. It wires together three active pages inside a
single `QStackedWidget` so the app feels like one continuous product rather
than separate screens. Uploading a dataset on the first page starts a model
run against the FastAPI backend (`backend/main.py`); the second page streams
its live progress and log output:

    1. `DatasetInitializationPage` — upload a dataset, configure a run.
    2. `ModelExecutionLogPage` — live pipeline status + log tail.
    3. `FinalResultsPage` — predicted drug rankings and patient profile.

`ModelVisualizationPage` (a decorative placeholder — fake canvas + hardcoded
run stats) is currently hidden from the navigation flow; see the `!Hidden`
comments below to restore it.

Every page shares the same visual language (sidebar, header, cards, tables)
via `styles.theme` and `widgets.navigation` / `widgets.tables`, so this file
is only responsible for constructing the pages and wiring navigation
callbacks between them — it contains no page-specific layout code.

How to run:
    1. Install dependencies (PySide6 at minimum): `pip install PySide6`.
       If a `requirements.txt` is present, prefer `pip install -r requirements.txt`.
    2. Run this file: `python main.py`.
    3. Assumes it is executed from the project root so the `pages/`,
       `widgets/`, and `styles/` packages are importable; no other project
       structure is assumed.
"""

from __future__ import annotations
import sys

from PySide6.QtWidgets import QApplication, QMainWindow, QStackedWidget

from pages.dataset_initialization_page import DatasetInitializationPage
from pages.final_results_page import FinalResultsPage
from pages.model_execution_log_page import ModelExecutionLogPage

# !Hidden: ModelVisualizationPage is a decorative placeholder (fake canvas +
# hardcoded run stats + a frozen progress bar), hidden from the navigation
# flow until it's wired to real data.
# from pages.model_visualization_page import ModelVisualizationPage
from styles.theme import apply_theme


def main() -> int:
    """Construct the application window, wire page navigation, and run the event loop.

    Builds all four pages up front, adds them to a `QStackedWidget`, and
    connects each page's navigation callbacks to `stack.setCurrentWidget`
    closures so pages never need direct references to one another.

    Returns:
        The Qt application's exit code, suitable for `sys.exit`/`SystemExit`.
    """
    app = QApplication(sys.argv)
    apply_theme(app)

    window = QMainWindow()
    window.setWindowTitle("MSC16 - Dataset Initialization")

    stack = QStackedWidget()

    # Set whenever a run starts (`show_log_page`), read by
    # `show_final_results_page` — a simple shared reference so the run id
    # and target cell line don't need to be threaded through every page
    # callback signature.
    current_run: dict[str, str | None] = {"run_id": None, "target_cell_line": None}

    # Each page is built with callbacks pointing at these closures rather
    # than at each other, so pages stay decoupled from one another and only
    # need to know about the navigation *action*, not the destination page.
    def show_upload_page() -> None:
        stack.setCurrentWidget(upload_page)

    def show_log_page(run_info: dict) -> None:
        # Only for the upload page's explicit call (with the run-start
        # response) — never connect this directly to a Qt signal, since
        # PySide would pass the `clicked(checked: bool)` argument through
        # as `run_info`.
        current_run["run_id"] = run_info["run_id"]
        current_run["target_cell_line"] = run_info.get("target_cell_line")
        log_page.start_watching(run_info)
        stack.setCurrentWidget(log_page)

    def show_log_page_only() -> None:
        stack.setCurrentWidget(log_page)

    # !Hidden: ModelVisualizationPage is unreferenced (see the import above).
    # def show_visualization_page() -> None:
    #     stack.setCurrentWidget(visualization_page)

    def show_final_results_page() -> None:
        if current_run["run_id"] is not None:
            final_results_page.load_results(current_run["run_id"])
        if current_run["target_cell_line"] is not None:
            final_results_page.set_sample_id(current_run["target_cell_line"])
        stack.setCurrentWidget(final_results_page)

    upload_page = DatasetInitializationPage(on_initialize_upload=show_log_page)
    log_page = ModelExecutionLogPage(
        on_upload_clicked=show_upload_page,
        # !Hidden: no on_model_visualization_clicked — see the import above.
        on_run_complete=show_final_results_page,
    )
    # !Hidden: ModelVisualizationPage is unreferenced (see the import above).
    # visualization_page = ModelVisualizationPage(
    #     on_upload_clicked=show_upload_page,
    #     on_model_running_clicked=show_log_page_only,
    #     on_model_analytics_clicked=show_log_page_only,
    #     on_finish_clicked=show_final_results_page,
    # )
    final_results_page = FinalResultsPage(
        on_upload_clicked=show_upload_page,
        on_model_running_clicked=show_log_page_only,
    )

    stack.addWidget(upload_page)
    stack.addWidget(log_page)
    # !Hidden: visualization_page is unreferenced (see the import above).
    # stack.addWidget(visualization_page)
    stack.addWidget(final_results_page)
    stack.setCurrentWidget(upload_page)

    window.setCentralWidget(stack)
    window.resize(1280, 1024)
    window.show()

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
