# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
#Go up one level to find the project
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cladeomatic'
copyright = '2023, James Robertson'
author = 'James Robertson'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#To install sphinx: pip install sphinx
#To install the read the docs theme: pip install sphinx-rtd-theme

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/logo.png'
html_title = "Clade-o-matic"
globaltoc_includehidden = True
toc_object_entries_show_parents = 'all'
add_module_names = True
html_sidebars = {
    '**': [
        'globaltoc.html',
        'relations.html',
        'searchbox.html'
    ]
}

html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'navigation_depth': -1,
}

########### TRICK FOUND ON SOME TUTORIAL : ADD IN THE MOCK_MODULES ANY EXTERNAL MODULE YOU'RE USING IN YOUR PACKAGE.

#import mock

#MOCK_MODULES = ['pandas', 'scipy', 'ray', 'matplotlib', 'Bio', 'deprecate']
#for mod_name in MOCK_MODULES:
#    sys.modules[mod_name] = mock.Mock()