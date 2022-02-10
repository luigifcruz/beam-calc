import pandas as pd

from logic import Beams
from io import StringIO
from bokeh.layouts import column, layout, row
from cartography import CoordEcef
from bokeh.tile_providers import get_provider, ESRI_IMAGERY
from bokeh.plotting import figure, curdoc, gridplot
from bokeh.models import LinearColorMapper, ColorBar, Button
from bokeh.models.widgets import DataTable, TableColumn, TextAreaInput, TextInput
from bokeh.events import PointEvent


def generate_data():
    global beams
    global phi_d_scalar
    global theta_d_scalar
    global frequency_scalar
    global antenna_map_src
    global coordinates_text
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

    antenna_map_src.data = dict(latlon_coords)

    beams = Beams(
        frequency_scalar.value,
        phi_d_scalar.value,
        theta_d_scalar.value,
        enu_coords)

    l = beams.azimuth[-1]
    r = beams.azimuth[0]
    b = beams.zenith[-1]
    t = beams.zenith[0]

    beam_pattern_plot_src.data = {
        "image": [beams.db.T],
        "x":     [r],
        "y":     [t],
        "dw":    [l-r],
        "dh":    [b-t],
    }

    update_powerplot_data(-1, -1)


def update_powerplot_data(azi, zen):
    global beams
    global zenith_plot_src
    global azimuth_plot_src

    if azi == -1 and azi == -1:
        azi = beams.db.shape[1] // 2
        zen = beams.db.shape[0] // 2

    azimuth_plot_src.data = {
        "angle": beams.azimuth - beams.theta_d,
        "power": beams.db[:, azi],
    }

    zenith_plot_src.data = {
        "angle": beams.zenith - beams.phi_d,
        "power": beams.db[zen],
    }


def beam_pattern_plot():
    fig = figure(
        title="Beam Pattern",
        x_axis_label="Azimuth Angle [deg]",
        y_axis_label="Zenith Angle [deg]",
        x_axis_location="above",
        tools="crosshair",
        sizing_mode='scale_both',
    )

    color_mapper = LinearColorMapper(palette='Viridis256', low=-30.0)

    plot = fig.image(color_mapper=color_mapper)

    def on_mouse_move(event: PointEvent):
        global beams
        global beam_pattern_plot_src

        y = beam_pattern_plot_src.data["y"][0]
        dw = beam_pattern_plot_src.data["dw"][0]
        mw = beams.db.shape[1] - 1

        azi = int(((event.y - y) / dw) * mw)
        azi = azi if azi < mw else mw        

        x = beam_pattern_plot_src.data["x"][0]
        dh = beam_pattern_plot_src.data["dh"][0]
        mh = beams.db.shape[0] - 1

        zen = int(((event.x - x) / dh) * mh)
        zen = zen if zen < mh else mh      

        update_powerplot_data(azi, zen)

    def on_mouse_leave():
        update_powerplot_data(-1, -1)

    fig.on_event('mouseleave', on_mouse_leave)
    fig.on_event('mousemove', on_mouse_move)

    return (fig, plot.data_source)


def zenith_plot():
    fig = figure(
        x_axis_label="Power [dB]",
        y_axis_label="Zenith Angle Offset [deg]",
        y_axis_location="right",
        width=225,
    )

    plot = fig.line(x="power", y="angle")

    return (fig, plot.data_source)


def azimuth_plot():
    fig = figure(
        y_axis_label="Power [dB]",
        x_axis_label="Azimuth Angle Offset [deg]",
        y_axis_location="right",
        height=200,
    )

    plot = fig.line(x="angle", y="power")

    return (fig, plot.data_source)


def plots_interface():
    global zenith_plot_src
    global azimuth_plot_src
    global beam_pattern_plot_src

    zenith_plot_obj, zenith_plot_src = zenith_plot()
    azimuth_plot_obj, azimuth_plot_src = azimuth_plot()
    beam_pattern_plot_obj, beam_pattern_plot_src = beam_pattern_plot()

    return gridplot([
        [beam_pattern_plot_obj, zenith_plot_obj],
        [azimuth_plot_obj, None]
    ])


def coordinates_interface():
    global coordinates_text
    global coordinates_table_src

    with open('ata-datum-condensed.csv', 'r') as file:
        default_ecef = file.read()

    coordinates_text = TextAreaInput(
        value=default_ecef,
        title="Array Coordinates (ECEF)",
        rows=25,
    )

    button_obj = Button(label="Update Coordinates")
    button_obj.on_click(generate_data)

    return column([
        coordinates_text, 
        button_obj,
    ], sizing_mode='scale_width')


def antenna_map():
    fig = figure(
        title="Array Antennas (WGS-84)",
        x_axis_label="Longitude [deg]",
        y_axis_label="Latitude [deg]",
        x_axis_type="mercator",
        y_axis_type="mercator",
        sizing_mode='scale_both',
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


def scalars_interface():
    global antenna_map_src
    global frequency_scalar
    global phi_d_scalar
    global theta_d_scalar

    def on_change(attr, old, new):
        del attr, old, new
        generate_data()

    frequency_scalar = TextInput(value="500", title="Frequency (MHz)")
    frequency_scalar.on_change("value", on_change)

    phi_d_scalar = TextInput(value="30", title="Zenith (deg)")
    phi_d_scalar.on_change("value", on_change)

    theta_d_scalar = TextInput(value="30", title="Azimuth (deg)")
    theta_d_scalar.on_change("value", on_change)

    antenna_map_obj, antenna_map_src = antenna_map()

    return column([
        frequency_scalar,
        row([
            phi_d_scalar,
            theta_d_scalar,
        ]),
        antenna_map_obj,
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
    ], sizing_mode='scale_width')

    doc = curdoc()
    doc.add_root(plots)

main()
