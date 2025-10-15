# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------

project = 'SPACER'
copyright = '2025, Jia Yao and Tao Wang'
author = 'Jia Yao'

release = '1.0'
version = '1.0.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex',
    'nbsphinx',
]

# Let Sphinx read Jupyter notebooks but not execute them
nbsphinx_allow_errors = True
nbsphinx_execute = 'never'

# Ignore some notebook warnings
nbsphinx_prolog = """
.. raw:: html

   <style>
       .nbinput.nblast, .nboutput.nblast { margin-bottom: 1em; }
   </style>
"""

# -- Intersphinx -------------------------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_logo = "_static/spacer_logo.png"
html_favicon = "_static/favicon.ico"

html_theme_options = {
    'logo_only': True,
    'display_version': False,
}

# -- Options for EPUB output -------------------------------------------------
epub_show_urls = 'footnote'
