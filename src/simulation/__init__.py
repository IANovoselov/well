from flask import Flask
app = Flask(__name__)
from .entrypoints import flask_app
