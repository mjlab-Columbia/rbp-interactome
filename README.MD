# Getting Started

This is an interactive tool for viewing RBP interaction data collected from IP and SEC mass spectrometry methods.

## Installation

Start off by cloning the repository localy. Open up terminal and navigate to the folder where the library is present and download the required libraries with: 

```bash
pip3 install -r requirements.txt
```

## Running the web server

After installing the necessary requirements run this command to start the web server:

```bash
bokeh serve web_app.py
```
You should see an output simillar to this:

![alt text](https://github.com/aa4793/rbp_analysis/blob/main/host.jpg?raw=true)

Copy and paste the url from this output into your browser in our case:

```bash
http://localhost:5006/web_app
```
and you should see the web appplication:

![alt text](https://github.com/aa4793/rbp_analysis/blob/main/web.jpg?raw=true)
