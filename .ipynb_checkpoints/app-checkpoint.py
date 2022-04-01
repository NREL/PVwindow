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

centeredstyle={'width':'10%','background-color':'red','position': 'absolute','top':'50%','-ms-transform':'translateY(-50%)','transform':'translateY(-50%)'}
listingstyle = {'width':'50%','float':'left','padding':'6px 0','font-size':'16px'}

#{'display':'inline-block',
# 'font-size':'22px','width':'15%',
# 'height':'36px',"verticalAlign": "center",
# 'text-align':'center','float':'left'}

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
            html.H2('Device Stack',className=''),#,class_name="me-md-1"

            html.Div(
                    [
                        dbc.Row(
                            dbc.Col(html.H4("Parameters:"))
                        ),
                        #First input row
                        dbc.Row(
                            [
                                dbc.Col(html.Div(["R",html.Sub('series'), ' [Ω]']),width={'size':3,'offset':1}),
                                dbc.Col(dbc.Input(id='Rs_in'),width={'size':2,'offset':0}),
                                dbc.Col(html.Div(["R",html.Sub('shunt'), ' [Ω]']),width={'size':3,'offset':0}),
                                dbc.Col(dbc.Input(id='Rsh_in'),width={'size':2,'offset':0}),
                            ],
                            justify='start',
                            align='center'
                        ),
                        #2nd input row
                        dbc.Row(
                            [
                                dbc.Col(html.Div(["U",html.Sub("outside"), ' [Wm',html.Sup('-2'),'K',html.Sup('-1'),']']),
                                        width={'size':3,'offset':1}),
                                dbc.Col(dbc.Input(id='Uout_in'),width={'size':2,'offset':0}),
                                dbc.Col(html.Div(["U",html.Sub("inside"), ' [Wm',html.Sup('-2'),'K',html.Sup('-1'),']']),
                                        width={'size':3,'offset':0}),
                                dbc.Col(dbc.Input(id='Uin_in'),width={'size':2,'offset':0}),
                            ],
                            justify='start',
                            align='center'
                        ),
                        #3rd input row
                        dbc.Row(
                            [
                                dbc.Col(html.Div(["T",html.Sub("outside"), ' [℃]']),width={'size':3,'offset':1}),
                                dbc.Col(dbc.Input(id='Tout_in'),width={'size':2,'offset':0}),
                                dbc.Col(html.Div(["T",html.Sub("inside"), ' [℃]']),width={'size':3,'offset':0}),
                                dbc.Col(dbc.Input(id='Tin_in'),width={'size':2,'offset':0}),
                            ],
                            justify='start',
                            align='center'
                        ),
                        #4th input row
                        dbc.Row(
                            [
                                dbc.Col(html.Div(["θ",html.Sub("inc"), ' [deg]']),width={'size':3,'offset':1}),
                                dbc.Col(dbc.Input(id='theta_in'),width={'size':2,'offset':0}),
                            ],
                            justify='start',
                            align='center'
                        ),
                        dbc.Row(
                            dbc.Col(html.H4("Layers:",className='')),
                        ),
                        dbc.Row(
                            [
                                
                                dbc.Col(dbc.Button('+',
                                                color="primary",
                                                id='button_add',
                                                n_clicks=0,
                                                className='col-12'),
                                        width={'size':2,'offset':1}
                                       ),
                                dbc.Col(dbc.ButtonGroup([
                                     dbc.Button('↑',
                                                outline=True,color="primary",
                                                id='button_up',n_clicks=0),
                                     dbc.Button('↓',
                                                outline=True,color="primary",
                                                id='button_down',n_clicks=0),
                                     ],
                                 className='col-12'
                                 ),width={'size':3,'offset':0})
                                
                            ],
                            align='center',#style={'font-size':'18px'},
                        ),
                        html.Br(),
                        dbc.Row(
                                [
                                    dbc.Col(dbc.Button('Solve',
                                                   id='button_solve',
                                                   color='secondary',
                                                   n_clicks=0,
                                                   className='col-12'),
                                           width={'size':6,'offset':2}
                                           ),
                                    dbc.Col(
                                                dbc.Spinner(html.Div(id='waiting-for-stuff',className='m-1'),color='primary'),

                                           )
                                ]
                        )
                        
                    ]
                ),

             html.Br(),
             #html.H1(''),
             dash_table.DataTable(
                 columns = [{'name':'Material','id':'Material',
                             'presentation': 'dropdown'},
                            {'name':'Thickness [μm]',
                             'id':'Thickness [μm]','presentation': 'input'},
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
                style_cell={'fontSize':16}, 
                            #'font-family':'Arial, Helvetica, sans-serif'},
                #workaround so bootstrap formatting works with dropdowns... see https://github.com/plotly/dash-table/issues/221
                css=[{"selector": ".Select-menu-outer", 
                      "rule": "display: block !important"}],
                 
             ),
            
        ],
        style = {'float':'left','width':'50%',
                 'border-right-style':'none','border-width':'thin'},
    ),
    #right side of interface
    dbc.Container(
        [
            html.H2('Analysis'),#,className='display-6'),
            dbc.Tabs(
            [
                dbc.Tab(label="A, R, T", tab_id="rat-tab"),#,
                        #tab_style={"marginLeft": "auto"}),
                dbc.Tab(label="Metrics", tab_id="metrics-tab"),
            ],
            id="tabs",
            active_tab="rat-tab",
            ),
            html.Div(id="tab-content",className="p-4"),
        ],
        style = {'float':'left','width':'50%','position':'relative'},
    ),
#    html.Div(
#        [
#            html.Div('first',style=centeredstyle)
#        ],
#        style={'height':'36px','width':'50%','float':'left','background-#color':'yellow','position':'relative'})
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
            VLT = data['metric-stuff']['VLT']
            PCE = data['metric-stuff']['PCE']
            SHGC = data['metric-stuff']['SHGC']
            Tcell = data['metric-stuff']['Tcell']
            
            thestuff = [
                dbc.Row(
                            [
                                dbc.Col(html.Div("Transmitted color:",),width={'size':5,'offset':0}),
                                dbc.Col(dbc.Button('   ',
                                                   className='p-3 col-12',
                                                   style={'background-color':color,'border-color':'black'}
                                                  ),
                                        width={'size':2,'offset':0}
                                       ),

                            ],
                            justify='start',
                            align='center'
                        ),
                dbc.Row(
                            [
                                dbc.Col(html.Div("Transmitted chromaticity:",),width={'size':5,'offset':0}),
                                dbc.Col(dbc.Button('   ',
                                                   className='p-3 col-12',
                                                   style={'background-color':chromaticity,'border-color':'black'}
                                                  ),
                                        width={'size':2,'offset':0}
                                       ),

                            ],
                            justify='start',
                            align='center'
                        ),
                dbc.Row(
                            [
                                dbc.Col(html.Div("Visible light transmission:",),width={'size':5,'offset':0}),
                                dbc.Col(dbc.Button(str(VLT)[:5],
                                                   className='p-1 col-12',
                                                   style={'color':'black','background-color':'transparent','border-color':'black'}
                                                  ),
                                        width={'size':2,'offset':0}
                                       ),

                            ],
                            justify='start',
                            align='center'
                        ),
                dbc.Row(
                            [
                                dbc.Col(html.Div("Power conversion efficiency:",),width={'size':5,'offset':0}),
                                dbc.Col(dbc.Button(str(PCE)[:5],
                                                   className='p-1 col-12',
                                                   style={'color':'black','background-color':'transparent','border-color':'black'}
                                                  ),
                                        width={'size':2,'offset':0}
                                       ),

                            ],
                            justify='start',
                            align='center'
                        ),
                dbc.Row(
                            [
                                dbc.Col(html.Div("Solar heat gain coefficient:",),width={'size':5,'offset':0}),
                                dbc.Col(dbc.Button(str(SHGC)[:5],
                                                   className='p-1 col-12',
                                                   style={'color':'black','background-color':'transparent','border-color':'black'}
                                                  ),
                                        width={'size':2,'offset':0}
                                       ),

                            ],
                            justify='start',
                            align='center'
                        ),
                dbc.Row(
                            [
                                dbc.Col(html.Div("Cell Temperature:",),width={'size':5,'offset':0}),
                                dbc.Col(dbc.Button([str(Tcell)[:5],' ℃'],
                                                   className='p-1 col-12',
                                                   style={'color':'black','background-color':'transparent','border-color':'black'}
                                                  ),
                                        width={'size':2,'offset':0}
                                       ),

                            ],
                            justify='start',
                            align='center'
                        ),
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
    Output('waiting-for-stuff','children'),
    Input('button_solve','n_clicks'),
    State('layer_table', 'data'),
    State('Rs_in','value'),
    State('Rsh_in','value'),
    State('Uin_in','value'),
    State('Uout_in','value'),
    State('Tin_in','value'),
    State('Tout_in','value'),
    State('theta_in','value'),
)
def make_figures(n_clicks,data_from_table,Rs,Rsh,Uin,Uout,Tin,Tout,theta):
        
    
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
            #family = "Open Sans",
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
    metrics = stack.get_transmitted_color(lambdas,i_ang)
    #print(metrics)
    
    ready = []
    if Rsh:
        print('Rsh from input: ' + Rsh)   
        ready.append(True)
    else:
        ready.append(False)
    if Rs:
        print('Rs from input: ' + Rs)    
        ready.append(True)
    else:
        ready.append(False)
        
    if Uin:
        print('Uin from input: ' + Uin)    
        ready.append(True)
    else:
        ready.append(False)
        
    if Uout:
        print('Uout from input: ' + Uout)
        ready.append(True)
    else:
        ready.append(False)
        
    if Tin:
        print('Tin from input: ' + Tin)
        ready.append(True)
    else:
        ready.append(False)
        
    if Tout:
        print('Tout from input: ' + Tout)
        print(type(Tout))
        print(type(float(Tout)))
        ready.append(True)
    else:
        ready.append(False)
        
    if theta:
        print('theta from input: ' + theta)
        ready.append(True)
    else:
        ready.append(False)
        
    
    if all(ready):
        print('let us go')
        thegoods = wpv.get_performance_characteristics(stack,float(Tin)+273.15,float(Tout)+273.15,float(Uin),float(Uout),
                                                        float(Rs),float(Rsh),float(theta))
        metrics['PCE']=thegoods['PCE']
        metrics['SHGC']=thegoods['SHGC']
        metrics['Tcell']=thegoods['Tcell']-273.15
        #thegoods = wpv.get_performance_characteristics_old(stack,1,float(Tin),float(Tout),float(Uin),float(Uout),
        #                                                float(Rs),float(Rsh),1,float(theta))
        print('did it')
        print(thegoods)
                                                    
    else:
        metrics['PCE']='0.0'
        metrics['SHGC']='N/A'
        metrics['Tcell']='N/A'
        print('more inputs needed')
    
    VLT = stack.get_visible_light_transmission(lambdas,i_ang)
    metrics['VLT']=VLT
    print(metrics)
                            
        
    return {'rat-fig':fig,'metric-stuff':metrics},None


if __name__ == '__main__':
    app.run_server(debug=True)