import pandas as pd
import numpy as np
import os
from time import sleep
import matplotlib.pyplot as plt
from plot import plot_CNV
from plot_helpers import convert_base


###### cnv files #######


def load_plot_df(sample, folder="", ext="cluster"):
    """
    load the combines SNP/cov data file for visualization
    """
    file = os.path.join(folder, sample, f"{sample}.{ext}")
    plot_df = pd.read_csv(file, sep="\t")
    return plot_df


def load_cnv_df(sample, folder, cnv_tool="cnacs", edit_suffix="edit", ascat_penalty=25):
    """
    loads the appropriate CNV detection tool output and creates a df from it
    also the tool-specific output_file is returned
    """

    # here comes the switch for the different tools
    cnv_folder = os.path.join(folder, sample)

    file = os.path.join(cnv_folder, f"{sample}.cnv.{cnv_tool}.txt")
    if cnv_tool in ["ascat", "refphase"]:
        file = file.replace(".txt", f"{ascat_penalty}.txt")

    file_edited = file.replace(".txt", f".{edit_suffix}.txt")

    # if edited file is existing, use it
    if os.path.isfile(file_edited):
        cnv_df = pd.read_csv(file_edited, sep="\t")
        return cnv_df, file_edited

    # else, use the original files and re-format
    cnv_df = pd.read_csv(file, sep="\t")
    cnv_df["Call"] = "not_validated"

    # adjust cols to CNACS standard
    if cnv_tool in ["ascat", "refphase"]:
        cnv_df = cnv_df.rename(
            dict(
                chrom="Chr",
                start="Start",
                end="End",
                sample_id="sample",
                cn_a="CNmajor",
                cn_major="CNmajor",
                cn_b="CNminor",
                cn_minor="CNminor",
            ),
            axis=1,
        )
        cnv_df["Chr"] = "chr" + cnv_df["Chr"].astype(str)
        cnv_df.loc[cnv_df["Chr"] == "chr23", "Chr"] = "ChrX"
        cnv_df["sample"] = sample
    cols = ["sample"] + [
        col for col in cnv_df.columns if not col in ["sample", "hg19", "liftoverInfo"]
    ]

    return cnv_df.loc[:, cols], file_edited


def load_cnv_data(
    sample,
    CNV_folder="",
    ascat_penalty=25,
    edit_suffix="edit",
    raw_extension="cluster",
    cnv_tool="cnacs",
):
    """
    setup a validation for a sample
    """

    plot_df = load_plot_df(sample, folder=CNV_folder, ext=raw_extension)

    cnv_df, out_file = load_cnv_df(
        sample,
        folder=CNV_folder,
        edit_suffix=edit_suffix,
        cnv_tool=cnv_tool,
        ascat_penalty=ascat_penalty,
    )

    return plot_df, cnv_df, out_file


##### VIZ update functions  ###############
def get_region(row, zoom=1, padding=0.25):
    """
    takes a region string in short form chr3:4Mb-6Mb and returns a
    chr3:40000000-60000000
    """
    start = row["Start"]
    end = row["End"]
    dist = end - start
    center = start + (dist / 2)
    shift = dist * (padding + zoom / 2)
    start, end = max(0, int(center - shift)), int(center + shift)
    return f"{row['Chr']}:{start}-{end}"


def get_row_output(row, tool="cnacs", view="all", zoom=1):
    """
    prepares for every row of cnacs_data and returns the output dict
    """
    out = {}
    out["chrom"] = row["Chr"] if view == "chrom" else "all"
    out["region"] = get_region(row, zoom=zoom) if view == "region" else ""
    out["marker_region"] = f"{row['Chr']}:{row['Start']}-{row['End']}"

    # get generic info field information
    start = f"{round(row['Start'] / 1e6, 1)}Mb"
    end = f"{round(row['End'] / 1e6, 1)}Mb"
    sample_info = f"Sample {row['sample']}@{row['Chr']}:{start}-{end} | "
    out["info"] = f"Call: {row['Call']} | {sample_info}"

    # get tool-specific info
    if tool == "cnacs":
        out["info"] += f"TCN = {round(row['TCN'],2)} | {row['SNPcount']} SNPs "
        if not np.isnan(row["BAF"]):
            out["info"] += f" | {round(row['BAF'], 2)}"
    # for ASCAT or refphase, only use the CN major/minor info
    else:
        out["info"] += f"CNmajor={row['CNmajor']} CNminor={row['CNminor']}"
    return out


def update_plot(plot_df, cnv_df, axes, counter, tool="cnacs", fig_params={}):
    for ax in axes:
        ax.cla()
    row = cnv_df.iloc[counter, :]
    # print(row)
    # get the tool specific output
    out = get_row_output(
        row, tool=tool, view=fig_params["view"], zoom=fig_params["zoom"]
    )

    if out["chrom"] != "all":
        out["chrom"] = [out["chrom"]]
    axes = plot_CNV(
        plot_df,
        axes,
        info=out["info"],
        chroms=out["chrom"],
        region=out["region"],
        marker_region=out["marker_region"],
        **fig_params,
    )


from ipywidgets import widgets, HBox, VBox


def make_handlers(
    plot_df,
    cnv_df,
    *,
    sample_file="",
    out_file="",
    tool="cnacs",
    axes,
    counter,
    fig_params,
):
    def fwd(b):
        nonlocal counter
        counter += 1
        if counter >= len(cnv_df.index):
            counter = 0
        update_plot(plot_df, cnv_df, axes, counter, tool=tool, fig_params=fig_params)

    def back(b):
        nonlocal counter
        if counter == 0:
            counter = len(cnv_df.index) - 1
        else:
            counter -= 1
        update_plot(plot_df, cnv_df, axes, counter, tool=tool, fig_params=fig_params)

    def all_view(b):
        if fig_params["view"] != "all":
            fig_params["view"] = "all"
            update_plot(
                plot_df, cnv_df, axes, counter, tool=tool, fig_params=fig_params
            )

    def chrom_view(b):
        if fig_params["view"] != "chrom":
            fig_params["view"] = "chrom"
            update_plot(
                plot_df, cnv_df, axes, counter, tool=tool, fig_params=fig_params
            )

    def region_view(b):
        if fig_params["view"] != "region":
            fig_params["view"] = "region"
            update_plot(
                plot_df, cnv_df, axes, counter, tool=tool, fig_params=fig_params
            )

    def call_ok(b):
        cnv_df.loc[counter, "Call"] = "OK"
        fwd(b)

    def call_reject(b):
        cnv_df.loc[counter, "Call"] = "Rejected"
        update_plot(plot_df, cnv_df, axes, counter, tool=tool, fig_params=fig_params)
        fwd(b)

    def save(b):
        cnv_df.to_csv(out_file, sep="\t", index=False)
        print(f"File saved to {out_file}")

    def exclude(b):
        # get sample name from cnv_df
        sample = cnv_df.loc[0, "sample"]
        sample_df = pd.read_csv(sample_file, sep="\t")
        sample_df.loc[sample_df["sample"] == sample, "exclude"] = True
        sample_df.to_csv(sample_file, sep="\t", index=False)
        print(f"Sample {sample} excluded")

    handlers = dict(
        fwd=fwd,
        back=back,
        all_view=all_view,
        chrom_view=chrom_view,
        region_view=region_view,
        call_ok=call_ok,
        call_reject=call_reject,
        save=save,
        exclude=exclude,
    )
    return handlers


def make_widget(
    sample,
    CNV_folder="",
    sample_file="",
    edit_suffix="edit",
    raw_extension="cluster",
    cnv_tool="cnacs",
    ascat_penalty=25,
    min_snp=4,
    min_range=0,
    fig_params={},
    button_desc=dict(
        back="Back",
        fwd="Forward",
        all_view="AllView",
        chrom_view="ChromView",
        region_view="RegionView",
        call_ok="AcceptCall",
        call_reject="RejectCall",
        save="SaveResults",
        exclude="ExcludeSample",
    ),
):

    # load the data
    plot_df, cnv_df, out_file = load_cnv_data(
        sample,
        CNV_folder=CNV_folder,
        edit_suffix=edit_suffix,
        raw_extension=raw_extension,
        cnv_tool=cnv_tool,
        ascat_penalty=ascat_penalty,
    )
    # combine the default fig params with additional params given as arguments
    fig_params_default = dict(
        colormap="coolwarm_r",
        color_chroms=True,
        ylims=dict(cov=(-1, 2.5), snp=(0, 1)),
        label_size=8,
        figsize=(10, 4),
        marker_alpha=0.3,
        zoom=1,
        view="chrom",
    )

    fig_params_default.update(fig_params)
    fig_params = fig_params_default.copy()

    # reduce the cnv_df using restrictions
    if cnv_tool == "cnacs":
        cnv_df = cnv_df.query("SNPcount >= @min_snp")
    else:
        cnv_df = cnv_df.query("CNmajor != 1 or CNminor != 1")

    min_range = convert_base(min_range)
    # use the min_range and reset to make counter work
    cnv_df = cnv_df.query("End - Start > @min_range").reset_index(drop=True)

    # prepare the figure and init the counter
    fig, axes = plt.subplots(
        2, figsize=fig_params["figsize"], gridspec_kw={"height_ratios": [1, 2]}
    )
    counter = 0

    # show plot for row0
    update_plot(plot_df, cnv_df, axes, counter, tool=cnv_tool, fig_params=fig_params)

    # create the handlers
    handlers = make_handlers(
        plot_df,
        cnv_df,
        sample_file=os.path.join(CNV_folder, sample_file),
        out_file=out_file,
        counter=counter,
        tool=cnv_tool,
        axes=axes,
        fig_params=fig_params,
    )

    ## create the buttons
    buttons = {}
    for button in button_desc:
        buttons[button] = widgets.Button(description=button_desc[button])

    # assign the handlers
    for button in buttons:
        buttons[button].on_click(handlers[button])

    # arrange the buttons
    button_list = list(buttons.values())
    nav_buttons = HBox(button_list[:2])
    view_buttons = HBox(button_list[2:5])
    validate_buttons = HBox(button_list[5:])
    return VBox([nav_buttons, view_buttons, validate_buttons])


###################################### apply CNV ################################


def retrieve_CNV_data(sample, CNV_path, save_cnv=True):
    """
    gathers the edited CNV data for one sample and returns the CNV_df and the info from the sample_list
    if all.txt already exists, it just loads it
    """

    # read sample info from sample file
    sample_info = (
        pd.read_csv(os.path.join(CNV_path, "CNV_sample_list2.txt"), sep="\t")
        .set_index("sample")
        .loc[sample, :]
    )

    # get base folder
    folder = os.path.join(CNV_path, sample)

    # looks for the all file and if found, returns it with the sample info
    all_file = os.path.join(folder, f"{sample}.cnv.all.txt")
    if os.path.isfile(all_file):
        print(f"Found combined CNV data file {all_file}.")
        cnv_df = pd.read_csv(all_file, sep="\t")
        return cnv_df, sample_info

    # find all edit files in the sample folder
    edit_files = [file for file in next(os.walk(folder))[2] if "edit" in file]
    # concat the cnv dfs
    dfs = []
    for file in edit_files:
        print(f"Found edit-file {file}.")
        df = pd.read_csv(os.path.join(folder, file), sep="\t")
        dfs.append(df)
        # get the tool name from the file name
        df["tool"] = file.split(".")[2]
    cnv_df = pd.concat(dfs)
    # sort by chrom and coord
    cnv_df["Chr"] = pd.Categorical(
        cnv_df["Chr"], [f"chr{chrom}" for chrom in list(range(1, 23)) + ["X", "Y"]]
    )
    cnv_df = cnv_df.sort_values(["Chr", "Start"]).reset_index(drop=True)

    # save the gathered file
    if save_cnv:
        print(f"Saving combined CNV data file to {all_file}")
        cnv_df.to_csv(all_file, sep="\t", index=False)
    return cnv_df, sample_info


def CNV2filter(df, row):
    """"""

    chrom, start, end = row.loc[["Chr", "Start", "End"]]
    # get the boolean mask for assigning
    cnv_muts = (df["Chr"] == chrom) & (df["Start"] >= start) & (df["End"] <= end)
    df.loc[cnv_muts, "CNV"] = True
    ####
    # here the other info can be applied as well


def apply_CNV(filter_df, *, sample, CNV_path):
    """
    takes a cohort_filter_df and applies the CNV_data to it
    """
    cnv_df, sample_info = retrieve_CNV_data(sample, CNV_path)
    # reduce cnv_df to OK calls
    ok_df = cnv_df.query('Call == "OK"')

    # reduce to sample
    f_df = filter_df.query("sample == @sample").copy()

    # apply the CNV info
    f_df.loc[:, "CNV"] = False
    for _, row in ok_df.iterrows():
        CNV2filter(f_df, row)

    return f_df, cnv_df