import numpy as np
import pandas as pd

from logic import Beams

from io import StringIO
from bokeh.layouts import column, layout, row
from cartography import CoordEcef
from bokeh.tile_providers import get_provider, ESRI_IMAGERY
from bokeh.plotting import figure, curdoc, gridplot
from bokeh.models import LinearColorMapper, ColorBar, Button
from bokeh.models.widgets import DataTable, TableColumn, TextAreaInput, TextInput


def generate_data():
    global frequency_scalar
    global phi_d_scalar
    global theta_d_scalar

    global antenna_map_src
    global coordinates_text
    global coordinates_table_src
    global zenith_plot_src
    global azimuth_plot_src
    global beam_pattern_plot_src

    ecef_coords = pd.read_csv(StringIO(str(coordinates_text.value)))
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
        frequency_scalar.value,
        phi_d_scalar.value,
        theta_d_scalar.value,
        enu_coords)

    azimuth_plot_src.data = {
        "angle": np.rad2deg(beams.azimuth - beams.theta_d),
        "power": beams.db[:, beams.db.shape[1]//2],
    }

    zenith_plot_src.data = {
        "angle": np.rad2deg(beams.zenith - beams.phi_d),
        "power": beams.db[beams.db.shape[0]//2],
    }

    l, r, b, t = np.rad2deg([
        beams.azimuth[-1],
        beams.azimuth[0],
        beams.zenith[-1],
        beams.zenith[0]
    ])

    beam_pattern_plot_src.data = {
        "image": [beams.db.T],
        "x": [r],
        "y": [b],
        "dw": [l-r],
        "dh": [b-t],
    }

    coordinates_table_src.data = dict(ecef_coords)
    antenna_map_src.data = dict(latlon_coords)


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


def beam_pattern_plot():
    fig = figure(
        title="Beam Pattern",
        x_axis_label="Azimuth [deg]",
        y_axis_label="Zenith Angle [deg]"
    )

    color_mapper = LinearColorMapper(palette='Viridis256', low=-30.0)

    plot = fig.image(color_mapper=color_mapper)

    fig.add_layout(
        ColorBar(
            color_mapper=color_mapper,
            title="Beam Gain [dB]"
        ),
        place='below'
    )

    return (fig, plot.data_source)


def zenith_plot():
    fig = figure(
        title=f"Zenith Angle Offset",
        x_axis_label="Zenith Angle Offset [deg]",
        y_axis_label="Power [dB]"
    )

    plot = fig.line(x="angle", y="power")

    return (fig, plot.data_source)


def azimuth_plot():
    fig = figure(
        title=f"Azimuth Angle Offset",
        x_axis_label="Azimuth Offset [deg]",
        y_axis_label="Power [dB]")

    plot = fig.line(x="angle", y="power")

    return (fig, plot.data_source)


def plots_interface():
    global antenna_map_src
    global zenith_plot_src
    global azimuth_plot_src
    global beam_pattern_plot_src

    antenna_map_obj, antenna_map_src = antenna_map()
    zenith_plot_obj, zenith_plot_src = zenith_plot()
    azimuth_plot_obj, azimuth_plot_src = azimuth_plot()
    beam_pattern_plot_obj, beam_pattern_plot_src = beam_pattern_plot()

    return gridplot([
        [antenna_map_obj, beam_pattern_plot_obj],
        [azimuth_plot_obj, zenith_plot_obj]
    ], sizing_mode="scale_height")


def coordinates_table():
    table = DataTable(columns=[
        TableColumn(title="Antenna", field="antenna"),
        TableColumn(title="X", field="x"),
        TableColumn(title="Y", field="y"),
        TableColumn(title="Z", field="z")
    ])

    return (table, table.source)


def coordinates_interface():
    global coordinates_text
    global coordinates_table_src

    with open('ata-datum-condensed.csv', 'r') as file:
        default_ecef = file.read()

    coordinates_text = TextAreaInput(
        value=default_ecef,
        title="Array Coordinates (ECEF)"
    )

    coordinates_table_obj, coordinates_table_src = coordinates_table()

    button_obj = Button(label="Update Coordinates")
    button_obj.on_click(generate_data)

    return column([
        row([
            coordinates_text, 
            coordinates_table_obj,
        ], sizing_mode="scale_width"),
        button_obj,
    ], sizing_mode="scale_width")


def scalars_interface():
    global frequency_scalar
    global phi_d_scalar
    global theta_d_scalar

    frequency_scalar = TextInput(value="500", title="Frequency (MHz)")
    phi_d_scalar = TextInput(value="30", title="Zenith (deg)")
    theta_d_scalar = TextInput(value="30", title="Azimuth (deg)")

    return column([
        frequency_scalar,
        phi_d_scalar,
        theta_d_scalar,
    ])


def main():
    global antenna_map_src

    coordinates_interface_obj = coordinates_interface()
    scalars_interface_obj = scalars_interface()
    plots_interface_obj = plots_interface()

    generate_data()

    plots = layout(children=[
        row([
            scalars_interface_obj,
            plots_interface_obj,
        ]),
        coordinates_interface_obj,
    ], sizing_mode="scale_width")
    curdoc().add_root(plots)

main()
