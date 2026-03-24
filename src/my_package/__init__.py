# src/my_package/__init__.py
# This file makes Python treat this directory as a package
from src.my_package.core import add

__all__ = ['add']