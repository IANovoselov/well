from flask import Flask, jsonify, request
from .. import app


@app.route("/index", methods=["POST", "GET"])
def index():
    gain = request.json["gain"]
    time_constant = request.json["gain"]
    if gain and time_constant:
        result = step_first_order_lag(gain, time_constant)
        return jsonify(result), 200
    return "OK", 200
