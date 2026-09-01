"""
Centralized path resolution for SpecHLA.

Priority chain:
  1. SPECHLA_DB env var (user override)
  2. $CONDA_PREFIX/share/spechla/db (conda install)
  3. No implicit database fallback

Usage:
  from spechla_paths import get_db_dir, get_script_dir
"""

__all__ = [
    "get_db_dir",
    "get_default_db_dir",
    "get_imgt_dir",
    "get_assets_dir",
    "get_script_dir",
]

import os

_dev_path_done = False

def _setup_dev_path():
    """In dev mode, add repo's bin/ to PATH so vendored binaries are found."""
    global _dev_path_done
    if _dev_path_done:
        return
    _dev_path_done = True
    conda = os.environ.get('CONDA_PREFIX', '')
    if conda and os.path.isdir(os.path.join(conda, 'share', 'spechla')):
        return
    bin_dir = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin'))
    if not os.path.isdir(bin_dir):
        return
    path_dirs = os.environ.get('PATH', '').split(os.pathsep)
    if bin_dir in path_dirs:
        return
    dirs = [bin_dir,
            os.path.join(bin_dir, 'SpecHap', 'build'),
            os.path.join(bin_dir, 'extractHairs', 'build'),
            os.path.join(bin_dir, 'fermikit', 'fermi.kit')]
    os.environ['PATH'] = os.pathsep.join(dirs) + os.pathsep + os.environ.get('PATH', '')

def get_db_dir():
    """Return the built reference directory, raising if none is available."""
    _setup_dev_path()
    env = os.environ.get('SPECHLA_DB', '')
    if env and os.path.isdir(os.path.join(env, 'HLA')):
        return env
    conda = os.environ.get('CONDA_PREFIX', '')
    conda_db = os.path.join(conda, 'share', 'spechla', 'db')
    if conda and os.path.isdir(os.path.join(conda_db, 'HLA')):
        return conda_db
    raise RuntimeError(
        'SpecHLA reference not found. Set SPECHLA_DB to a directory built by '
        'script/build_reference.py.'
    )


def get_default_db_dir():
    """Return where a reference build should be written by default.

    Unlike get_db_dir() this does not require the directory to exist yet, so
    the reference builder can use it as the default output location.
    """
    env = os.environ.get('SPECHLA_DB', '')
    if env:
        return env
    conda = os.environ.get('CONDA_PREFIX', '')
    if conda:
        return os.path.join(conda, 'share', 'spechla', 'db')
    return ''


def get_imgt_dir():
    """Return a local IMGT/HLA checkout, raising if none is configured."""
    env = os.environ.get('SPECHLA_IMGT', '')
    if env and os.path.isdir(env):
        return env
    conda = os.environ.get('CONDA_PREFIX', '')
    conda_imgt = os.path.join(conda, 'share', 'spechla', 'imgt')
    if conda and os.path.isdir(conda_imgt):
        return conda_imgt
    raise RuntimeError(
        'IMGT/HLA source not found. Set SPECHLA_IMGT to a local IMGT/HLA clone.'
    )


def get_assets_dir():
    """Return the bundled reference construction assets shipped with SpecHLA."""
    env = os.environ.get('SPECHLA_ASSETS', '')
    if env and os.path.isdir(env):
        return env
    conda = os.environ.get('CONDA_PREFIX', '')
    conda_assets = os.path.join(conda, 'share', 'spechla', 'reference_assets')
    if conda and os.path.isdir(conda_assets):
        return conda_assets
    repo_assets = os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'share', 'reference_assets'))
    if os.path.isdir(repo_assets):
        return repo_assets
    raise RuntimeError(
        'Reference construction assets not found. Set SPECHLA_ASSETS to the '
        'directory holding extend.fa and representative_alleles.tsv.'
    )


def get_script_dir():
    """Return the directory holding the SpecHLA helper scripts."""
    _setup_dev_path()
    env = os.environ.get('SPECHLA_SCRIPT', '')
    if env and os.path.isdir(env):
        return env
    conda = os.environ.get('CONDA_PREFIX', '')
    conda_script = os.path.join(conda, 'share', 'spechla', 'script')
    if conda and os.path.isdir(conda_script):
        return conda_script
    return os.path.dirname(os.path.abspath(__file__))
