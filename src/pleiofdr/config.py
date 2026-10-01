"""Text configuration files (port of TextConfig.m and the defaults declared in runme.m)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

# Defaults exactly as declared in runme.m. Values are stored as strings, like TextConfig.declare
# (which calls num2str); booleans become '1'/'0'.
DEFAULTS: dict[str, str] = {
    "traitfolder": "../example_data_for_pleiotropy",
    "traitfile1": "PGC2_SCZ.mat",
    "traitname1": "SCZ",
    "traitfiles": "{'COG_charge.mat'}",
    "traitnames": "{'COGNITION'}",
    "reffile": "ref9545380_1kgPhase3eur_LDr2p1.mat",
    "refinfo": "",
    "mlibrary": "./",
    "randprune": "1",
    "randprune_gc": "0",
    "reset_pruneidx": "1",
    "randprune_n": "20",
    "randprune_file": "",
    "randprune_repeats": "default",
    "stattype": "conjfdr",
    "fdrthresh": "0.05",
    "pthresh": "1",
    "onscreen": "0",
    "outputdir": "test",
    "manh_fontsize_genenames": "12",
    "manh_legend": "NorthEast",
    "manh_plot": "1",
    "manh_yspace": "0.75",
    "manh_ymargin": "0.25",
    "manh_colorlist": "[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]",
    "exclude_chr_pos": "[6 25119106 33854733]",
    "exclude_from_discovery": "0",
    "mafthresh": "0.005",
    "use_standard_gc": "0",
    "perform_gc": "1",
    "exclude_ambiguous_snps": "0",
    "dummy_zscore": "0",
    "exit_matlab_upon_completion": "0",
}

_SCALARS = {
    "true": 1.0,
    "false": 0.0,
    "nan": np.nan,
    "inf": np.inf,
    "+inf": np.inf,
    "-inf": -np.inf,
    "-nan": np.nan,
}


class ConfigError(ValueError):
    pass


def _parse_scalar(token: str, key: str) -> float:
    lowered = token.lower()
    if lowered in _SCALARS:
        return _SCALARS[lowered]
    try:
        return float(token)
    except ValueError:
        raise ConfigError(f'config "{key}" has an invalid number: "{token}"') from None


def parse_matrix(text: str, key: str = "") -> np.ndarray:
    """Parse a MATLAB numeric literal such as '0.05', '[6 1 2; 8 3 4]' or 'true' (str2num subset).

    Returns a 2-D float array; an empty string or '[]' gives shape (0, 0).
    """
    body = text.strip()
    if body.startswith("["):
        if not body.endswith("]"):
            raise ConfigError(f'config "{key}" is not valid: "{text}"')
        body = body[1:-1]
    rows = [r for r in re.split(r"[;\n]", body)]
    parsed = [[_parse_scalar(t, key) for t in re.split(r"[\s,]+", r.strip()) if t] for r in rows]
    parsed = [r for r in parsed if r]
    if not parsed:
        return np.zeros((0, 0))
    if len({len(r) for r in parsed}) != 1:
        raise ConfigError(f'config "{key}" is not valid: "{text}" (ragged rows)')
    return np.array(parsed, dtype=float)


def parse_cell(text: str, key: str = "") -> list[str]:
    """Parse a MATLAB cell array of strings such as "{'a.mat', 'b.mat'}"."""
    body = text.strip()
    if not (body.startswith("{") and body.endswith("}")):
        raise ConfigError(f'config "{key}" is not a cell')
    body = body[1:-1]
    items: list[str] = []
    i = 0
    while i < len(body):
        ch = body[i]
        if ch in " \t,;":
            i += 1
            continue
        if ch not in "'\"":
            raise ConfigError(f'config "{key}" is not a cell of strings: "{text}"')
        quote, i, chars = ch, i + 1, []
        while True:
            if i >= len(body):
                raise ConfigError(f'config "{key}" has an unterminated string: "{text}"')
            if body[i] == quote:
                if i + 1 < len(body) and body[i + 1] == quote:  # '' escapes a quote
                    chars.append(quote)
                    i += 2
                    continue
                i += 1
                break
            chars.append(body[i])
            i += 1
        items.append("".join(chars))
    return items


class TextConfig:
    """Key/value configuration; keys must be declared before a file can set them."""

    def __init__(self, defaults: dict[str, str] | None = None) -> None:
        self.values: dict[str, str] = dict(DEFAULTS if defaults is None else defaults)

    def declare(self, key: str, default: str) -> None:
        self.values[key] = default

    def load_file(self, filename: str | Path) -> None:
        print(f'reading config file "{filename}"')
        for raw in Path(filename).read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("=")
            if len(parts) < 2:
                raise ConfigError(f'invalid config line (expected key=value): "{line}"')
            # TextConfig.m keeps only the text between the first and second '='
            key, value = parts[0].strip(), parts[1].strip()
            print(f"loaded {key}={value}")
            self._check(key)
            self.values[key] = value

    def _check(self, key: str) -> None:
        if key not in self.values:
            raise ConfigError(f'config has no member "{key}"')

    def get_str(self, key: str) -> str:
        self._check(key)
        return self.values[key]

    def get_num(self, key: str) -> float:
        v = parse_matrix(self.get_str(key), key)
        if v.shape != (1, 1):
            raise ConfigError(f'config "{key}" has size > 1')
        return float(v[0, 0])

    def get_bool(self, key: str) -> bool:
        v = parse_matrix(self.get_str(key), key)
        # MATLAB: ~(v==0); an empty value is false in an if-statement
        return bool(v.size) and bool(np.all(v != 0))

    def get_mat(self, key: str) -> np.ndarray:
        return parse_matrix(self.get_str(key), key)

    def get_cell(self, key: str) -> list[str]:
        return parse_cell(self.get_str(key), key)

    def print(self) -> None:
        print("configuration:")
        for key, value in self.values.items():
            print(f"{key}\t=\t{value}")
