import numpy as np
import pandas as pd

from logic import Beams

from io import StringIO
from bokeh.layouts import column, layout, row
from cartography import CoordEcef
from bokeh.tile_providers import ESRI_IMAGERY, get_provider
from bokeh.plotting import figure, curdoc, gridplot
from bokeh.models import LinearColorMapper, ColorBar, ColumnDataSource, Button
from bokeh.models.widgets import DataTable, TableColumn, TextAreaInput, TextInput, Slider


def generate_data():
    global frequency_scalar_src
    global phi_d_scalar_src
    global theta_d_scalar_src

    global antenna_map_src
    global coordinates_text_src
    global coordinates_table_src

    ecef_coords = pd.read_csv(StringIO(str(coordinates_text_src.value)))
    latlon_coords = pd.DataFrame(columns=['antenna', 'lat', 'lon', 'alt'])
    enu_coords = pd.DataFrame(columns=['antenna', 'e', 'n', 'u'])

    ref_ant = ecef_coords.iloc[0]
    ref_utm = CoordEcef(ref_ant.x, ref_ant.y, ref_ant.z).to_utm()

    for _, row in ecef_coords.iterrows():
        ant_utm = CoordEcef(row.x, row.y, row.z).to_utm()

        pos = ant_utm.to_latlon().to_webmercator()
        coord = pd.DataFrame({
            "antenna": row.antenna,
            "lat": pos.lat,
            "lon": pos.lon,
            "alt": pos.alt,
        }, index=[0])
        latlon_coords = pd.concat([latlon_coords, coord])

        pos = ant_utm - ref_utm
        coord = pd.DataFrame({
            "antenna": row.antenna,
            "e": pos.e,
            "n": pos.n,
            "u": pos.u,
        }, index=[0])
        enu_coords = pd.concat([enu_coords, coord])

    beams = Beams(
        frequency_scalar_src.value,
        phi_d_scalar_src.value,
        theta_d_scalar_src.value,
        enu_coords)

    coordinates_table_src.data = dict(ColumnDataSource(data=ecef_coords).data)
    antenna_map_src.data = dict(ColumnDataSource(data=latlon_coords).data)
    print("updated")


def antenna_map():
    fig = figure(
        title="Array Antennas (WGS-84)",
        x_axis_label="Longitude [deg]",
        y_axis_label="Latitude [deg]",
        x_axis_type="mercator",
        y_axis_type="mercator"
    )

    color_mapper = LinearColorMapper(palette='Viridis256')

    plot = fig.circle(
        x="lat",
        y="lon",
        size=10.0,
        color={"field": "alt", "transform": color_mapper}
    )

    fig.add_layout(
        ColorBar(
            color_mapper=color_mapper,
            title="Altitude [m]"
        ),
        place='below'
    )

    tile_provider = get_provider(ESRI_IMAGERY)
    fig.add_tile(tile_provider)

    return (fig, plot.data_source)


def coordinates_text():
    with open('ata-datum-condensed.csv', 'r') as file:
        default_ecef = file.read()

    text = TextAreaInput(value=default_ecef, title="ECEF CSV File")

    return (text, text)


def coordinates_table():
    table = DataTable(columns=[
        TableColumn(title="Antenna", field="antenna"),
        TableColumn(title="X", field="x"),
        TableColumn(title="Y", field="y"),
        TableColumn(title="Z", field="z")
    ])

    return (table, table.source)


def coordinates_interface():
    global coordinates_table_src
    global coordinates_text_src

    coordinates_text_obj, coordinates_text_src = coordinates_text()
    coordinates_table_obj, coordinates_table_src = coordinates_table()

    button_obj = Button(label="Update Coordinates")
    button_obj.on_click(generate_data)

    return column([
        row([
            coordinates_text_obj, 
            coordinates_table_obj,
        ], sizing_mode="scale_width"),
        button_obj,
    ], sizing_mode="scale_width")


def scalars_interface():
    global frequency_scalar_src
    global phi_d_scalar_src
    global theta_d_scalar_src

    frequency_scalar_src = TextInput(value="500", title="Frequency (MHz)")
    phi_d_scalar_src = Slider(start=0.0, end=90.0, value=30.0, title="Zenith (deg)")
    theta_d_scalar_src = Slider(start=0.0, end=360.0, value=30.0, title="Azimuth (deg)")

    return column([
        frequency_scalar_src,
        phi_d_scalar_src,
        theta_d_scalar_src,
    ])

def main():
    global antenna_map_src

    antenna_map_obj, antenna_map_src = antenna_map()
    coordinates_interface_obj = coordinates_interface()
    scalars_interface_obj = scalars_interface()

    generate_data()

    plots = layout(children=[
        row([
            scalars_interface_obj,
            antenna_map_obj,
        ]),
        coordinates_interface_obj,
    ], sizing_mode="scale_width")
    curdoc().add_root(plots)

main()
