from __future__ import annotations

ICON_MAP: dict[str, str] = {
    "biotech": "⚗",
    "add": "+",
    "settings": "⚙",
    "history": "🕘",
    "help_outline": "?",
    "notifications": "🔔",
    "account_circle": "◉",
    "keyboard_arrow_down": "▾",
    "cloud_upload": "⇪",
    "search": "⌕",
    "insights": "◍",
    "terminal": "⌨",
    "download": "⤓",
    "hub": "◎",
    "analytics": "◔",
    "folder_data": "🗁",
    "architecture": "▦",
    "stop_circle": "⏹",
    "settings_system_daydream": "◉",
    "support": "❓",
    "summarize": "🧾",
    "info": "ⓘ",
    "open_in_new": "↗",
    "scatter_plot": "∷",
    "description": "🗎",
    "sync": "↻",
}


def icon_text(name: str) -> str:
    return ICON_MAP.get(name, name)
