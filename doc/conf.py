"""rm -rf doc/_generated/; python setup.py build_sphinx -E -a
"""

project = "easistrain"
release = "0.1"
copyright = "2021, ESRF"
author = "ESRF"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]
templates_path = ["_templates"]
exclude_patterns = []

always_document_param_types = True

html_theme = "classic"
html_static_path = []

autosummary_generate = True
autodoc_default_flags = [
    "members",
    "undoc-members",
    "show-inheritance",
]
