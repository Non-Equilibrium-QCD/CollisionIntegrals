# Configuration file for the Sphinx documentation builder.

project = 'CollIntegral'
copyright = '2025'
author = ''

# -- General configuration ---------------------------------------------------

extensions = [
    'myst_parser',
    'sphinx.ext.mathjax',
    'breathe',
]

# Breathe configuration
breathe_projects = {
    "CollIntegral": "../xml"
}
breathe_default_project = "CollIntegral"
breathe_default_members = ('members', 'undoc-members')

# Suppress duplicate C++ declaration warnings
suppress_warnings = ['duplicate_declaration.cpp']

# MyST configuration
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "colon_fence",
]

# Allow both .md and .rst
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# The master toctree document
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'furo'
html_static_path = ['_static']

# Furo theme options
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#2962ff",
        "color-brand-content": "#2962ff",
    },
    "dark_css_variables": {
        "color-brand-primary": "#82b1ff",
        "color-brand-content": "#82b1ff",
    },
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
}

# -- Options for LaTeX/PDF output --------------------------------------------

latex_engine = 'lualatex'

latex_elements = {
    'preamble': r'''
\usepackage{amsmath}
\usepackage{amssymb}
''',
    # Use Unicode fonts with LuaLaTeX
    'fontpkg': r'''
\usepackage{fontspec}
\setmonofont{DejaVu Sans Mono}[Scale=0.9]
''',
}

latex_documents = [
    (master_doc, 'CollIntegral.tex', 'CollIntegral Documentation',
     author, 'manual'),
]

# MathJax configuration for equation numbering
mathjax3_config = {
    'tex': {
        'tags': 'ams',
        'packages': {'[+]': ['ams']},
    },
}

# Number figures, tables, code-blocks
numfig = True
