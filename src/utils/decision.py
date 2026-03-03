from __future__ import annotations

def print_level(msg: str, level: str) -> None:
    level = level.upper()
    print(f"{level}: {msg}")

def ask_stop_if_abnormal(msg: str) -> bool:
    print(f"ABNORMAL: {msg}")
    ans = input("Stop now? (y/N) ").strip().lower()
    return ans == "y"