from dash import Dash, html, dcc, Input, Output
from dash import dash_table
#from dash.dependencies import Input, Output
#import dash_html_components as html
#import dash_core_components as dcc
import dash_bootstrap_components as dbc
#from dash.dash_table.Format import Group
#import dash_table
import pandas as pd
#import plotly.graph_objects as go
#import plotly.express as px

#print([{"name": i, "id": i} for i in df.columns])

app = Dash(external_stylesheets=[dbc.themes.SIMPLEX])
#app = Dash(__name__)

app.layout = html.Div(
    [
    html.Div(
        children=[
             html.Label('Window Configuration'),
             dash_table.DataTable(
                 data = [],
                 columns = [{'name':'Material','id':'Material'},
                            {'name':'Thickness [nm]','id':'Thickness [nm]'}],
             ),
             dbc.Button('+ Layer',id='button_add',outline=True,class_name='gap-4'),
             dbc.Button('↑ Layer',id='button_up'),
             dbc.Button('↓ Layer',id='button_down'),
             dbc.Button('- Layer',id='button_subtract'),
            
        ],
        style = {'float':'left','width':'40%'},
    ),
    
    html.Div(
        children=[
            html.Label('Tabs')
        ],
        style = {'float':'left','width':'60%'},
    )
    
    ]
 )


if __name__ == "__main__":
    app.run_server(debug=True)