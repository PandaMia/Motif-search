import json
import re
from pathlib import Path


def slugify_preset_name(name: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
    return slug or "preset"


def load_presets(base_dir: Path):
    presets = []
    data_dir = Path(base_dir, "data")

    for preset_path in sorted(data_dir.glob("*.json")):
        preset_data = json.loads(preset_path.read_text())
        stem = preset_path.stem
        dna = preset_data.get("DNA") or []

        presets.append(
            {
                "slug": preset_data.get("slug") or slugify_preset_name(stem),
                "name": preset_data.get("name") or stem,
                "filename": preset_path.name,
                "description": preset_data.get("description") or f"Preset loaded from {preset_path.name}.",
                "dna_count": len(dna) if isinstance(dna, list) else 0,
                "path": str(preset_path),
            }
        )

    return presets


def get_preset_map(base_dir: Path):
    return {preset["slug"]: preset for preset in load_presets(base_dir)}
