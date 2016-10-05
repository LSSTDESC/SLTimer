import sys
import os
import sphinx_rtd_theme

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',]

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

master_doc = 'index'
autosummary_generate = True
autoclass_content = "class"
autodoc_default_flags = ["members", "no-special-members"]

project = u'greeter'
author = u'Phil Marshall and the LSST Dark Energy Science Collaboration'
copyright = u'2016, ' + author
version = "0.1"
release = "0.1.0"
