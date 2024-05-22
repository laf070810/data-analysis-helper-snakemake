import ROOT

ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.lhcbStyle()
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT(snakemake.threads)
ROOT.gStyle.SetMarkerSize(0.5)


def draw_1Dhists(hists, output_file_path_prefix, output_file_path_suffix):
    for hist in hists:
        canvas = ROOT.TCanvas()
        hist.SetStats(False)
        hist.SetMarkerStyle(24)
        hist.SetMarkerColor(ROOT.kBlue)
        hist.SetLineColor(ROOT.kBlue)
        hist.SetMinimum(0.0)
        hist.SetXTitle(hist.GetTitle())
        hist.Draw("e")
        canvas.SaveAs(
            output_file_path_prefix + hist.GetName() + output_file_path_suffix
        )


def draw_2Dhists(hists, output_file_path_prefix, output_file_path_suffix):
    for hist in hists:
        canvas = ROOT.TCanvas()
        ROOT.gPad.SetLeftMargin(0.16)
        ROOT.gPad.SetRightMargin(0.16)
        hist.SetStats(False)
        hist.GetYaxis().SetTitleOffset(0.9)
        hist.Draw("colz")
        canvas.SaveAs(
            output_file_path_prefix + hist.GetName() + output_file_path_suffix
        )


for plot_info in snakemake.params.hists:
    datasets = plot_info["datasets"]
    tree_name = plot_info["tree_name"]
    hists_1D = plot_info.get("hists_1D", [])
    hists_2D = plot_info.get("hists_2D", [])
    df = ROOT.RDataFrame(tree_name, datasets)
    draw_1Dhists(
        [
            df.Filter(hist_1D.get("filter", "true"))
            .Define("__temp_var", hist_1D["expression"])
            .Histo1D(hist_1D["model"], "__temp_var")
            for hist_1D in hists_1D
        ],
        snakemake.params.output_file_path_prefix,
        snakemake.params.output_file_path_suffix,
    )
    draw_2Dhists(
        [
            df.Filter(hist_2D.get("filter", "true"))
            .Define("__temp_var_x", hist_2D["expression_x"])
            .Define("__temp_var_y", hist_2D["expression_y"])
            .Histo2D(hist_2D["model"], "__temp_var_x", "__temp_var_y")
            for hist_2D in hists_2D
        ],
        snakemake.params.output_file_path_prefix,
        snakemake.params.output_file_path_suffix,
    )
