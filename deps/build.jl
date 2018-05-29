# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("https://github.com/OptimalDesignLab/PkgFix.jl.git")
end
Pkg.checkout("PkgFix", "upgrade_0.6")

using PkgFix  # from now on, use PkgFix instead of Pkg for everything

pkgs = ["https://github.com/OptimalDesignLab/PumiInterface.jl.git" "upgrade_0.6";
#        "https://github.com/JuliaParallel/MPI.jl.git"  "v0.5.0";
        "https://github.com/OptimalDesignLab/SummationByParts.jl.git"  "jc_update_0.6";
#        "https://github.com/JuliaLang/ArrayViews.jl.git" "93e80390aeedb1dbcd90281b6dff7f760f430bc8";
#        "https://github.com/jipolanco/WriteVTK.jl.git"  "v0.6.1";
        "https://github.com/OptimalDesignLab/ODLCommonTools.jl.git" "upgrade_0.6"]



pkg_dict = PkgFix.installed()

for i=1:size(pkgs, 1)
  pkg_name = PkgFix.getRepoName(pkgs[i, 1])
  if !haskey(pkg_dict, pkg_name)
    PkgFix.add(pkgs[i, 1], branch_ish=pkgs[i, 2])
  end
end

global const ARRAYVIEWS_URL = "http://github.com/JuliaArrays/ArrayViews.jl.git"
global const ARRAYVIEWS_VER = "master"

if !haskey(pkg_dict, "ArrayViews")
  PkgFix.add(ARRYAVIEWS_URL, branch_ish=ARRAYVIEWS_VER)
else
  PkgFix.checkout("ArrayViews", ARRAYVIEWS_VER)
end



