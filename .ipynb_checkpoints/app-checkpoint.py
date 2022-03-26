#from jupyter_dash import JupyterDash
from dash import Dash
from dash import html, dcc, Input, Output, State
from dash import dash_table, callback_context
import numpy as np
#from dash.dependencies import Input, Output
#import dash_html_components as html
#import dash_core_components as dcc
import dash_bootstrap_components as dbc
#from dash.dash_table.Format import Group
#import dash_table
import pandas as pd
import plotly.graph_objects as go
import wpv
import os
#import plotly.express as px





matmap = {'Low Fe Glass':'nkLowFeGlass.csv',
          'FTO':'nkFTO.csv',
          'MAPI':'nkMAPI.csv',
          'Silver':'nkAg.csv',
          'AZO':'nkAZO.csv',
          'Bleached MAPI':'nkBleach.csv',
          'C60':'nkC60.csv',
          'ClAlPc':'nkClAlPc.csv',
          'EVA':'nkEVA.csv',
          'FTO':'nkFTO.csv',
          'ITO':'nkITO.csv',
          'MAPbBr₃':'nkMAPbBr3.csv',
          'NiO':'nkNiO.csv',
          'PTB7-Th:IEICO-4F':'nkPTB7_ThIEICO_4F.csv',
          'SiO₂':'nkSiO2.csv',
          'SnO₂ 1':'nkSnO2.csv',
          'SnO₂ 2':'nkSnO2_wBaseline.csv',
          'TiO₂':'nkTiO2.csv',
         }
matoptions = [{'label': i, 'value': matmap[i]} for i in matmap]

app = Dash(__name__,external_stylesheets=[dbc.themes.JOURNAL])
#app = JupyterDash(__name__,external_stylesheets=[dbc.themes.JOURNAL])
#app = JupyterDash(__name__,external_stylesheets=[dbc.themes.SANDSTONE])

#fixes issue with non initialized selection in table
app.config.suppress_callback_exceptions = True

# Create server variable with Flask server object for use with gunicorn
server = app.server

#ART_fig = 
Glass = wpv.Layer(4000,matmap['Low Fe Glass'],'i')
#FTO = wpv.Layer(0.3,'nkFTO','c')
#MAPI = wpv.Layer(0.06,'nkMAPI','c')
#Ag = wpv.Layer(0.01,'nkAg','c')
#TiO2lowE = wpv.Layer(0.02,'nkTiO2','c')
#EVA = wpv.Layer(1500,'nkEVA','i')


#layers = [Glass,FTO,MAPI,EVA,Glass,TiO2lowE,Ag,TiO2lowE]
layers = [Glass]

stack = wpv.Stack(layers)

num_lams = 200
lambdas = np.linspace(0.3,2.5,num=num_lams)
i_ang = 0
    
Is = stack.Is(lambdas)/np.max(stack.Is(lambdas))
LEF = stack.cieplf(lambdas)


app.layout = dbc.Container([
    dcc.Store(id="store"),
    dbc.Container([
        html.H4('PVwindow',className='display-4'),
        html.P(
                "A tool for analyzing semi-transparent photovoltaic devices",
                className="",
        ),
        html.Hr(className="my-2"),
        #html.Div('this thing drags',style={'width':'300px','height':'100px','background-color':'yellow'},draggable='true'),
        #html.Div('this thing drags',style={'width':'300px','height':'100px','background-color':'red'},draggable='true'),
    ]),
    # left side of interface  
    dbc.Container(
        [
             #html.H1('I am hungry'),
             html.H6('Device Stack',className='display-6'),#,class_name="me-md-1"
             html.Span('Layer:',style={'display':'inline-block','font-size':'22px','width':'20%','height':'36px',"verticalAlign": "center",'text-align':'center','float':'left'}),
                 dbc.Button('+',color="primary",id='button_add',style={'width':'10%','height':'36px','float':'left'},n_clicks=0),
             dbc.ButtonGroup([
                 dbc.Button('↑',outline=True,color="primary",id='button_up',n_clicks=0),
                 dbc.Button('↓',outline=True,color="primary",id='button_down',n_clicks=0),
             ],
             style={'width':'20%','height':'36px','float':'left'}
             ),
             #html.Br(),
             html.Div(style={'width':'5%','height':'36px','float':'left'}),
             dbc.Button('Solve',id='button_solve',style={'width':'20%','height':'36px'},color='secondary',n_clicks=0),#style={'width':'60%'},class_name='d-grid col-4 mx-auto'
             html.Br(),
             dash_table.DataTable(
                 columns = [{'name':'Material','id':'Material','presentation': 'dropdown'},
                            {'name':'Thickness [μm]','id':'Thickness [μm]','presentation': 'input'},
                            {'name':'PV','id':'PV','presentation': 'dropdown'}],
                 data=[{'Material':'nkLowFeGlass.csv',
                       'Thickness [μm]':4000,
                       'PV':False}],
                 id = 'layer_table',
                 editable=True,
                 row_selectable='single',
                 row_deletable=True,
                 selected_rows=[0],
                 style_as_list_view=True,
                 style_header={'text-align':'left'},
                 style_data={'text-align':'left'},
                 dropdown={
                    'Material': {
                        'options': matoptions
                    },
                    'PV': {
                         'options': [
                            {'label': 'Yes', 'value': True},
                            {'label': 'No', 'value': False}
                        ]
                    }
                },
                style_cell={'fontSize':16, 'font-family':'Arial, Helvetica, sans-serif'},
                #workaround so bootstrap formatting works with dropdowns... see https://github.com/plotly/dash-table/issues/221
                css=[{"selector": ".Select-menu-outer", 
                      "rule": "display: block !important"}],
                 
             ),
            
        ],
        style = {'float':'left','width':'50%','border-right-style':'none','border-width':'thin'},
    ),
    #right side of interface
    dbc.Container(
        [
            html.H6('Analysis',className='display-6'),
            dbc.Tabs(
            [
                dbc.Tab(label="A, R, T", tab_id="rat-tab",tab_style={"marginLeft": "auto"}),
                dbc.Tab(label="Metrics", tab_id="metrics-tab"),
            ],
            id="tabs",
            active_tab="rat-tab",
            ),
            html.Div(id="tab-content",className="p-4"),
        ],
        style = {'float':'left','width':'50%'},
    ),
    ],
    #fluid=True,
)


'''
The following leaned on this example: https://dash-bootstrap-components.opensource.faculty.ai/examples/graphs-in-tabs/
'''
@app.callback(
    Output("tab-content", "children"),
    [Input("tabs", "active_tab"), Input("store", "data")],
)
def render_tab_content(active_tab, data):
    """
    This callback takes the 'active_tab' property as input, as well as the
    stored graphs, and renders the tab content depending on what the value of
    'active_tab' is.
    """
    if active_tab and data is not None:
        if active_tab == "rat-tab":
            return dcc.Graph(figure=data["rat-fig"])
        
        elif active_tab == "metrics-tab":
            color = data['metric-stuff']['color']
            chromaticity = data['metric-stuff']['chromaticity']
            
            thestuff = [
                html.Span('Transmitted color: ',style={'display':'inline-block','font-size':'16px','width':'40%','height':'36px',"verticalAlign": "center",'text-align':'center','float':'left'}),
                dbc.Button('  ',style={'width':'100px','height':'36px','background-color':color,'border-color':'black'}),
                html.Br(),
                html.Span('Transmitted chromaticity: ',style={'display':'inline-block','font-size':'16px','width':'40%','height':'36px',"verticalAlign": "center",'text-align':'center','float':'left'}),
                dbc.Button('  ',style={'width':'100px','height':'36px','background-color':chromaticity,'border-color':'black'})
            ]
            
                            
            
            return thestuff
            
    return "No data to display yet"

@app.callback(
    Output('layer_table', 'data'),
    Output('layer_table', 'selected_rows'),
    Input('button_add', 'n_clicks'),
    Input('button_up','n_clicks'),
    Input('button_down','n_clicks'),
    State('layer_table', 'data'),
    State('layer_table', 'columns'),
    State('layer_table', 'dropdown'),
    State('layer_table', 'selected_rows')
)
def manipulate_stack(n_clicks_add, n_clicks_up, n_clicks_down, data, columns, ddown, srows):
    ctx = callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    #print('button id: ' + str(button_id))
    #display(ddown)
    #x = {c['id']: '' for c in columns}
    #print('x=' + str(x))
    
    new_selection=srows
    chosen=srows[0]
    
    nlayer = np.shape(data)[0]
    #print(nlayer)
    
    if n_clicks_add > 0:
        if button_id == 'button_add':
            data.append({c['id']: '' for c in columns})
    
    if n_clicks_up > 0:
        if button_id == 'button_up':
            if nlayer>0 and chosen>0:
                #print('old srows: ' + str(srows))
                # popping both the elements from list
                first_ele = data.pop(chosen)  
                # inserting in each others positions
                data.insert(chosen-1, first_ele) 
                
                new_selection = [chosen-1]
                #print('new srows: ' + str(new_selection))
            else:
                print('nothing to move')
                
    if n_clicks_down > 0:
        if button_id == 'button_down':
            if nlayer>0 and chosen<nlayer:
                # popping both the elements from list
                first_ele = data.pop(chosen)  
                # inserting in each others positions
                data.insert(chosen+1, first_ele) 
                new_selection = [chosen+1]
            else:
                print('nothing to move')
    #print(data)
        
    
    #stack.self_summary()
    #print(new_selection)
    return data, new_selection

"""
put callback here whose input is the solve button (remove above) and whose output is a figure in the id='store' object... 
should return a dictionary with keys for each figure (rat-fig)
...
Now adding outputs for analysis tab: color, chromaticity, PCE, SHGC...
"""

@app.callback(
    Output('store','data'),
    Input('button_solve','n_clicks'),
    State('layer_table', 'data'),
)
def make_figures(n_clicks,data_from_table):
    
    fig = go.Figure()
    
    fig.add_trace(
        go.Scatter(
            x=lambdas,
            y=LEF,
            name='LEF',
            showlegend=True,
            line = dict(color='gray',width=1, dash='solid'),
        )
    )
    
    fig.add_trace(
        go.Scatter(
            x=lambdas,
            y=Is,
            name='Solar',
            showlegend=True,
            line = dict(color='goldenrod',width=1, dash='solid'),
        )
    )
    
    stack.update_from_dash(data_from_table)
    
    pvabs = stack.get_specular_PV_abs(lambdas,i_ang)
    fig.add_trace(
        go.Scatter(
            x=lambdas,
            y=pvabs,
            name='A<sub>PV</sub>',
            showlegend=True,
            line = dict(color='red',width=1, dash='dot'),
        )
    )

    spectra = stack.get_specular_RAT(lambdas,i_ang)
    spec_names = ['Reflectivity','Absorptivity','Transmissivity']
    colors = ['black','red','blue']
    for spectrum,spec_name,color in zip(spectra,spec_names,colors):
        fig.add_trace(
            go.Scatter(
                x = lambdas,
                y = spectrum,
                name=spec_name,
                showlegend=True,
                line = dict(color=color,width=3, dash='solid'),
            )
        )
    
        
    if n_clicks == 0:
        fig.for_each_trace(lambda trace: trace.update(visible="legendonly"))
    else:
        fig.for_each_trace(lambda trace: trace.update(visible=True))
        
    fig.update_layout(
        xaxis_title="Wavelength, um",
        margin=dict(l=0,r=0,t=90,b=0),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        font=dict(
            family = "Arial, Helvetica, sans-serif",
            size=14,
        ),
        plot_bgcolor = 'rgba(0,0,0,0)',
        paper_bgcolor = 'rgba(0,0,0,0)',
        legend_font_size = 12,
    )
    fig.update_xaxes(showline=True,
                 linewidth=0.5,
                 linecolor='black',
                 #mirror=True,
                 ticks='outside')

    fig.update_yaxes(showline=True,
                 linewidth=0.5,
                 linecolor='black',
                 #mirror=True,
                 ticks='outside')

    #stack.self_summary()
    colorstuff = stack.get_transmitted_color(lambdas,i_ang)
    #print(colorstuff)
        
    return {'rat-fig':fig,'metric-stuff':colorstuff}


if __name__ == '__main__':
    app.run_server(debug=True)