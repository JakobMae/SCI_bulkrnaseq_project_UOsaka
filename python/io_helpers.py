"""Small shared path helpers for the gene-X SCI analysis scripts."""

from pathlib import Path
import yaml


def load_paths(config_file: str = "config/paths.yml") -> dict:
    config_path = Path(config_file)
    if not config_path.exists():
        raise FileNotFoundError(
            f"Missing {config_file}. Copy config/paths_template.yml to "
            "config/paths.yml and edit local paths."
        )

    with config_path.open("r") as handle:
        return yaml.safe_load(handle)


def project_path(relative_path: str, paths: dict | None = None) -> Path:
    if paths is None:
        paths = load_paths()

    path = Path(relative_path)
    if path.is_absolute():
        return path

    return Path(paths["project"]["root"]) / path


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def require_file(path: str | Path, label: str) -> Path:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {label}\nExpected: {path}")
    return path
