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

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/solar.csv')

print(df.to_dict('records'))
moop
#print([{"name": i, "id": i} for i in df.columns])

#app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
app = Dash(__name__)

app.layout = html.Div(
    [
    html.Div(
        children=[
             html.Label('Window Configuration'),
             dash_table.DataTable(
                 data = [],
                 columns = [{'name':'Material','id':'Material'},
                            {'name':'Thickness [nm]','id':'Thickness [nm]'}],
             )
        ],
        style = {'float':'left','width':'40%','background-color': '#a1edcc'},
    ),
    
    html.Div(
        children=[
            html.Label('Column 2')
        ],
        style = {'float':'left','width':'60%','background-color': '#f497f1'},
    )
    
    ]
 )


if __name__ == "__main__":
    app.run_server(debug=True)