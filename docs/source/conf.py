# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'QUPS'
copyright = '2024, Thurston Brevett'
author = 'Thurston Brevett'
release = 'v1.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinxcontrib.matlab', 'sphinx.ext.autodoc', 'sphinx.ext.autosummary']
primary_domain = "mat"
matlab_src_dir = "/home/thurston/sandbox/qups"
autoclass_content = "class" # show top of class help text
autodoc_member_order = "bysource"
# matlab_auto_link = "all" # not working :(
matlab_auto_link = "basic"
matlab_short_links = True

# defaults
templates_path = ['_templates']
exclude_patterns = []

autodoc_default_options = {
    'members': True,
    'show-inheritance': True
}



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
