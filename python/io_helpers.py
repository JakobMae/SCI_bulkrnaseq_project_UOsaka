"""Shared input/output helper functions for the Gpnmb SCI analysis repo."""

from pathlib import Path
import yaml


def load_paths(config_file: str = "config/paths.yml") -> dict:
    config_path = Path(config_file)
    if not config_path.exists():
        raise FileNotFoundError(
            f"Missing config file: {config_file}. "
            "Copy config/paths_template.yml to config/paths.yml and edit local paths."
        )
    with config_path.open("r") as handle:
        return yaml.safe_load(handle)


def project_path(*parts: str, paths: dict | None = None) -> Path:
    if paths is None:
        paths = load_paths()
    return Path(paths["project"]["root"]).joinpath(*parts)


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def check_file_exists(path: str | Path, label: str | None = None) -> Path:
    path = Path(path)
    if not path.exists():
        name = label or str(path)
        raise FileNotFoundError(f"Missing required file: {name}\nExpected path: {path}")
    return path
