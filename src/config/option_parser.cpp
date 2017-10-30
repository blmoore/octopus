
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "option_parser.hpp"

#include <vector>
#include <array>
#include <regex>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "utils/path_utils.hpp"
#include "utils/memory_footprint.hpp"
#include "utils/string_utils.hpp"
#include "basics/phred.hpp"
#include "exceptions/user_error.hpp"
#include "config.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace octopus { namespace options {

fs::path resolve_path(const fs::path& path, const OptionMap& options);
void parse_config_file(const fs::path& config_file, OptionMap& vm, const po::options_description& options);

// boost::option cannot handle option dependencies so we must do our own checks
void conflicting_options(const OptionMap& vm, const std::string& opt1, const std::string& opt2);
void option_dependency(const OptionMap& vm, const std::string& for_what, const std::string& required_option);
void check_positive(const std::string& option, const OptionMap& vm);
void check_reads_present(const OptionMap& vm);
void check_region_files_consistent(const OptionMap& vm);
void check_trio_consistent(const OptionMap& vm);
void validate_caller(const OptionMap& vm);
void validate(const OptionMap& vm);

po::parsed_options run(po::command_line_parser& parser);


OptionMap parse_options(const int argc, const char** argv)
{
    po::options_description general("General");
    general.add_options()
    ("help,h", "Produce help message")
    
    ("version", "Output the version number")
    
    ("config",
     po::value<fs::path>(),
     "A config file, used to populate command line options")
    
    ("debug",
     po::value<fs::path>()->implicit_value("octopus_debug.log"),
     "Writes verbose debug information to debug.log in the working directory")
    
    ("trace",
     po::value<fs::path>()->implicit_value("octopus_trace.log"),
     "Writes very verbose debug information to trace.log in the working directory")

    ("fast",
     po::bool_switch()->default_value(false),
     "Turns off some features to improve runtime, at the cost of decreased calling accuracy."
     " Equivalent to '-a off -l minimal x 50`")

    ("very-fast",
     po::bool_switch()->default_value(false),
     "The same as fast but also disables inactive flank scoring")
    ;
    
    po::options_description backend("Backend");
    backend.add_options()
    ("working-directory,w",
     po::value<fs::path>(),
     "Sets the working directory")
    
    ("threads",
     po::value<int>()->implicit_value(0),
     "Maximum number of threads to be used, enabling this option with no argument lets the application"
     " decide the number of threads ands enables specific algorithm parallelisation")
    
    ("max-reference-cache-footprint,X",
     po::value<MemoryFootprint>()->default_value(*parse_footprint("500MB"), "500MB"),
     "Maximum memory footprint for cached reference sequence")
    
    ("target-read-buffer-footprint,B",
     po::value<MemoryFootprint>()->default_value(*parse_footprint("6GB"), "6GB"),
     "None binding request to limit the memory footprint of buffered read data")
    
    ("max-open-read-files",
     po::value<int>()->default_value(250),
     "Limits the number of read files that can be open simultaneously")
    ;
    
    po::options_description input("I/O");
    input.add_options()
    ("reference,R",
     po::value<std::string>()->required(),
     "FASTA format reference genome file to be analysed. Target regions"
     " will be extracted from the reference index if not provded explicitly")
    
    ("reads,I",
     po::value<std::vector<fs::path>>()->multitoken(),
     "Space-separated list of BAM/CRAM files to be analysed."
     " May be specified multiple times")
    
    ("reads-file,i",
     po::value<fs::path>(),
     "File containing a list of BAM/CRAM files, one per line, to be analysed")
    
    ("one-based-indexing",
     po::bool_switch()->default_value(false),
     "Notifies that input regions are given using one based indexing rather than zero based")
    
    ("regions,T",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-separated list of regions (chrom:begin-end) to be analysed."
     " May be specified multiple times")
    
    ("regions-file,t",
     po::value<fs::path>(),
     "File containing a list of regions (chrom:begin-end), one per line, to be analysed")
    
    ("skip-regions,K",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-separated list of regions (chrom:begin-end) to skip"
     " May be specified multiple times")
    
    ("skip-regions-file,k",
     po::value<fs::path>(),
     "File of regions (chrom:begin-end), one per line, to skip")
    
    ("samples,S",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-separated list of sample names to analyse")
    
    ("samples-file,s",
     po::value<fs::path>(),
     "File of sample names to analyse, one per line, which must be a subset of the samples"
     " that appear in the read files")
    
    ("ignore-unmapped-contigs",
     po::bool_switch()->default_value(false),
     "Ignore any contigs that are not present in the read files")
    
    ("pedigree",
     po::value<fs::path>(),
     "PED file containing sample pedigree")
    
    ("output,o",
     po::value<fs::path>(),
     "File to where output is written. If unspecified, calls are written to stdout")
    
    ("contig-output-order",
     po::value<ContigOutputOrder>()->default_value(ContigOutputOrder::asInReferenceIndex),
     "The order contigs should be written to the output")
    
    ("legacy",
     po::bool_switch()->default_value(false),
     "Outputs a legacy version of the final callset in addition to the native version")
    
    ("regenotype",
     po::value<fs::path>(),
     "VCF file specifying calls to regenotype, only sites in this files will appear in the"
     " final output")
    ;
    
    po::options_description transforms("Read transformations");
    transforms.add_options()
    ("read-transforms",
     po::value<bool>()->default_value(true),
     "Enable all read transformations")

    ("mask-low-quality-tails",
     po::value<int>()->implicit_value(3),
     "Masks read tail bases with base quality less than this")
    
    ("soft-clip-masking",
     po::value<bool>()->default_value(true),
     "Turn on or off soft clip base recalibration")

    ("soft-clip-mask-threshold",
     po::value<int>()->implicit_value(3),
     "Only soft clipped bases with quality less than this will be recalibrated, rather than all bases")
    
    ("mask-soft-clipped-boundary-bases",
     po::value<int>()->default_value(2),
     "Masks this number of adjacent non soft clipped bases when soft clipped bases are present")
    
    ("adapter-masking",
     po::value<bool>()->default_value(true),
     "Enable adapter detection and masking")
    
    ("overlap-masking",
     po::value<bool>()->default_value(true),
     "Enable read segment overlap masking")
    ;
    
    po::options_description filters("Read filtering");
    filters.add_options()
    ("read-filtering",
     po::value<bool>()->default_value(true),
     "Enable all read filters")
    
    ("consider-unmapped-reads",
     po::bool_switch()->default_value(false),
     "Allows reads marked as unmapped to be used for calling")
    
    ("min-mapping-quality",
     po::value<int>()->default_value(20),
     "Minimum read mapping quality required to consider a read for calling")
    
    ("good-base-quality",
     po::value<int>()->default_value(20),
     "Base quality threshold used by min-good-bases and min-good-base-fraction filters")
    
    ("min-good-base-fraction",
     po::value<double>()->implicit_value(0.5),
     "Base quality threshold used by min-good-bases filter")
    
    ("min-good-bases",
     po::value<int>()->default_value(20),
     "Minimum number of bases with quality min-base-quality before read is considered")
    
    ("allow-qc-fails",
     po::bool_switch()->default_value(false),
     "Filters reads marked as QC failed")
    
    ("min-read-length",
     po::value<int>(),
     "Filters reads shorter than this")
    
    ("max-read-length",
     po::value<int>(),
     "Filter reads longer than this")
    
    ("allow-marked-duplicates",
     po::bool_switch()->default_value(false),
     "Allows reads marked as duplicate in alignment record")
    
    ("allow-octopus-duplicates",
     po::bool_switch()->default_value(false),
     "Allows reads considered duplicates by octopus")
    
    ("no-secondary-alignments",
     po::bool_switch()->default_value(false),
     "Filters reads marked as secondary alignments")
    
    ("no-supplementary-alignments",
     po::bool_switch()->default_value(false),
     "Filters reads marked as supplementary alignments")
    
    ("consider-reads-with-unmapped-segments",
     po::bool_switch()->default_value(false),
     "Allows reads with unmapped template segmenets to be used for calling")
    
    ("consider-reads-with-distant-segments",
     po::bool_switch()->default_value(false),
     "Allows reads with template segmenets that are on different contigs")
    
    ("no-adapter-contaminated-reads",
     po::bool_switch()->default_value(false),
     "Filter reads with possible adapter contamination")
    
    ("disable-downsampling",
     po::bool_switch()->default_value(false),
     "Disables downsampling")
    
    ("downsample-above",
     po::value<int>()->default_value(1000),
     "Downsample reads in regions where coverage is over this")
    
    ("downsample-target",
     po::value<int>()->default_value(500),
     "The target coverage for the downsampler")
    ;
    
    po::options_description variant_generation("Candidate variant generation");
    variant_generation.add_options()
    ("raw-cigar-candidate-generator,g",
     po::value<bool>()->default_value(true),
     "Enable candidate generation from raw read alignments (CIGAR strings)")
    
    ("assembly-candidate-generator,a",
     po::value<bool>()->default_value(true),
     "Enable candidate generation using local re-assembly")
    
    ("source-candidates,c",
     po::value<std::vector<fs::path>>()->multitoken(),
     "Variant file paths containing known variants. These variants will automatically become candidates")
    
    ("min-base-quality",
     po::value<int>()->default_value(20),
     "Only bases with quality above this value are considered for candidate generation")
    
    ("min-supporting-reads",
     po::value<int>()->implicit_value(2),
     "Minimum number of reads that must support a variant if it is to be considered a candidate."
     " By default octopus will automatically determine this value")
    
    ("max-variant-size",
     po::value<int>()->default_value(2000),
     "Maximum candidate variant size to consider (in region space)")
    
    ("kmer-sizes",
     po::value<std::vector<int>>()->multitoken()
     ->default_value(std::vector<int> {10, 15, 20}, "10 15 20")->composing(),
     "Kmer sizes to use for local assembly")
    
    ("num-fallback-kmers",
     po::value<int>()->default_value(7),
     "How many local assembly fallback kmer sizes to use if the default sizes fail")
    
    ("fallback-kmer-gap",
     po::value<int>()->default_value(10),
     "The gap size used to generate local assembly fallback kmers")
    
    ("max-region-to-assemble",
     po::value<int>()->default_value(400),
     "The maximum region size that can be used for local assembly")
    
    ("max-assemble-region-overlap",
     po::value<int>()->default_value(200),
     "The maximum number of bases allowed to overlap assembly regions")
    
    ("force-assemble",
     po::bool_switch()->default_value(false),
     "Forces all regions to be assembled")
    
    ("assembler-mask-base-quality",
     po::value<int>()->default_value(10),
     "Aligned bases with quality less than this will be converted to reference before"
     " being inserted into the De Bruijn graph")
    
    ("min-kmer-prune",
     po::value<int>()->default_value(2),
     "Minimum number of read observations to keep a kmer in the assembly graph before bubble extraction")

    ("max-bubbles",
     po::value<int>()->default_value(30),
     "Maximum number of bubbles to extract from the assembly graph")
    
    ("min-bubble-score",
     po::value<double>()->default_value(2.0),
     "Minimum bubble score that will be extracted from the assembly graph")
    ;
    
    po::options_description haplotype_generation("Haplotype generation");
    haplotype_generation.add_options()
    ("max-haplotypes,x",
     po::value<int>()->default_value(200),
     "Maximum number of candidate haplotypes the caller may consider. If a region contains"
     " more candidate haplotypes than this then filtering is applied")
    
    ("haplotype-holdout-threshold",
     po::value<int>()->default_value(2500),
     "Forces the haplotype generator to temporarily hold out some alleles if the number"
     " of haplotypes in a region exceeds this threshold")
    
    ("haplotype-overflow",
     po::value<int>()->default_value(200000),
     "Regions with more haplotypes than this will be skipped")
    
    ("max-holdout-depth",
     po::value<int>()->default_value(20),
     "Maximum number of holdout attempts the haplotype generator can make before the region"
     " is skipped")
    
    ("extension-level",
     po::value<ExtensionLevel>()->default_value(ExtensionLevel::normal),
     "Level of haplotype extension. Possible values are: conservative, normal, optimistic, aggressive")
    ;
    
    po::options_description caller("Caller (general)");
    caller.add_options()
    ("caller,C",
     po::value<std::string>()->default_value("population"),
     "Which of the octopus callers to use")
    
    ("organism-ploidy,P",
     po::value<int>()->default_value(2),
     "All contigs with unspecified ploidies are assumed the organism ploidy")
    
    ("contig-ploidies,p",
     po::value<std::vector<ContigPloidy>>()->multitoken()
     ->default_value(std::vector<ContigPloidy> {
        {boost::none, "Y", 1}, {boost::none, "MT", 1}}, "Y=1 MT=1")
     ->composing(),
     "Space-separated list of contig (contig=ploidy) or sample contig"
     " (sample:contig=ploidy) ploidies")
    
    ("contig-ploidies-file",
     po::value<fs::path>(),
     "File containing a list of contig (contig=ploidy) or sample contig"
     " (sample:contig=ploidy) ploidies, one per line")
    
    ("min-variant-posterior",
     po::value<Phred<double>>()->default_value(Phred<double> {2.0}),
     "Report variant alleles with posterior probability (phred scale) greater than this")
    
//    ("min-refcall-posterior",
//     po::value<Phred<double>>()->default_value(Phred<double> {2.0}),
//     "Report reference alleles with posterior probability (phred scale) greater than this")
//
//    ("refcall,v",
//     po::value<RefCallType>()->implicit_value(RefCallType::blocked),
//     "Caller will report reference confidence calls for each position (positional),"
//     " or in automatically sized blocks (blocked)")
    
    ("sites-only",
     po::bool_switch()->default_value(false),
     "Only reports call sites (i.e. without sample genotype information)")
    
    ("snp-heterozygosity,z",
     po::value<float>()->default_value(0.001, "0.001"),
     "SNP heterozygosity for the given samples")

    ("snp-heterozygosity-stdev",
     po::value<float>()->default_value(0.01, "0.01"),
     "Standard deviation of the SNP heterozygosity used for the given samples")
    
    ("indel-heterozygosity,y",
     po::value<float>()->default_value(0.0001, "0.0001"),
     "Indel heterozygosity for the given samples")
    
    ("use-uniform-genotype-priors",
    po::bool_switch()->default_value(false),
    "Use a uniform prior model when calculating genotype posteriors")

    ("explain-read-directions",
     po::value<bool>()->default_value(true),
     "Model read directions in the genotype likelihood calculation")
    ;
    
    po::options_description cancer("Caller (cancer)");
    cancer.add_options()
    ("normal-sample,N",
     po::value<std::string>(),
     "Normal sample - all other samples are considered tumour")
    
    ("somatic-mutation-rate",
     po::value<float>()->default_value(1e-05, "1e-05"),
     "Expected somatic mutation rate, per megabase pair, for this sample")
    
    ("min-expected-somatic-frequency",
     po::value<float>()->default_value(0.05, "0.05"),
     "Minimum expected somatic allele frequency in the sample")
    
    ("min-credible-somatic-frequency",
     po::value<float>()->default_value(0.01, "0.01"),
     "Minimum credible somatic allele frequency that will be reported")
    
    ("credible-mass",
     po::value<float>()->default_value(0.99, "0.99"),
     "Mass of the posterior density to use for evaluating allele frequencies")
    
    ("min-somatic-posterior",
     po::value<Phred<double>>()->default_value(Phred<double> {0.5}),
     "Minimum posterior probability (phred scale) to emit a somatic mutation call")
    
    ("somatics-only",
     po::bool_switch()->default_value(false),
     "Only report somatic variant calls")
    ;
    
    po::options_description trio("Caller (trio)");
    trio.add_options()
    ("maternal-sample,M",
     po::value<std::string>(),
     "Maternal sample")
    
    ("paternal-sample,F",
     po::value<std::string>(),
     "Paternal sample")
    
    ("snv-denovo-mutation-rate",
     po::value<float>()->default_value(1e-9, "1e-9"),
     "SNV de novo mutation rate, per base per generation")
    
    ("indel-denovo-mutation-rate",
     po::value<float>()->default_value(1e-10, "1e-10"),
     "INDEL de novo mutation rate, per base per generation")
    
    ("min-denovo-posterior",
     po::value<Phred<double>>()->default_value(Phred<double> {0.5}),
     "Minimum posterior probability (phred scale) to emit a de novo mutation call")
    
    ("denovos-only,d",
     po::bool_switch()->default_value(false),
     "Only report de novo variant calls (i.e. alleles unique to the child)")
    ;
    
    po::options_description phasing("Phasing");
    phasing.add_options()
    ("phasing-level,l",
     po::value<PhasingLevel>()->default_value(PhasingLevel::normal),
     "Level of phasing - longer range phasing can improve calling accuracy at the cost"
     " of runtime speed. Possible values are: minimal, conservative, moderate, normal, aggressive")
    
    ("min-phase-score",
     po::value<Phred<double>>()->default_value(Phred<double> {20.0}),
     "Minimum phase score (phred scale) required to report sites as phased")
    
    ("use-unconditional-phase-score",
     po::bool_switch()->default_value(false),
     "Computes unconditional phase scores rather than conditioning on called genotypes")
    ;
    
    po::options_description advanced("Advanced calling algorithm");
    advanced.add_options()
    ("haplotype-extension-threshold,e",
     po::value<Phred<double>>()->default_value(Phred<double> {100.0}, "100"),
     "Haplotypes with posterior probability less than this can be filtered before extension")
    
    ("inactive-flank-scoring",
     po::value<bool>()->default_value(true),
     "Disables additional calculation to adjust alignment score when there are inactive"
     " candidates in haplotype flanking regions")
    
    ("model-mapping-quality",
     po::value<bool>()->default_value(true),
     "Include the read mapping quality in the haplotype likelihood calculation")
    
    ("max-joint-genotypes",
     po::value<int>()->default_value(1000000),
     "The maximum number of joint genotype vectors to consider when computing joint"
     " genotype posterior probabilities")
    
    ("sequence-error-model",
     po::value<std::string>()->default_value("hiseq"),
     "The sequencer error model to use.")
    ;
    
    po::options_description call_filtering("Callset filtering");
    call_filtering.add_options()
    ("call-filtering",
     po::value<bool>()->default_value(false),
     "Enable all variant call filtering")
    
    ("model-filtering",
     po::value<bool>(),
     "Enable model based filtering of variant calls")
    ;
    
    po::options_description all("octopus options");
    all.add(general).add(backend).add(input).add(transforms).add(filters)
    .add(variant_generation).add(haplotype_generation).add(caller)
    .add(advanced).add(cancer).add(trio).add(phasing).add(call_filtering);
    
    OptionMap vm_init;
    po::store(run(po::command_line_parser(argc, argv).options(general).allow_unregistered()), vm_init);
    
    if (vm_init.count("help") == 1) {
        po::store(run(po::command_line_parser(argc, argv).options(caller).allow_unregistered()), vm_init);
        if (vm_init.count("caller") == 1) {
            const auto caller = vm_init.at("caller").as<std::string>();
            validate_caller(vm_init);
            if (caller == "individual") {
                std::cout << all << std::endl;
            } else if (caller == "population") {
                std::cout << all << std::endl;
            } else if (caller == "cancer") {
                std::cout << all << std::endl;
            } else {
                std::cout << all << std::endl;
            }
        } else {
            std::cout << all << std::endl;
        }
        return vm_init;
    }
    
    if (vm_init.count("version") == 1) {
        std::cout << "octopus " << config::Version << std::endl;
        return vm_init;
    }
    
    OptionMap vm;
    
    if (vm_init.count("config") == 1) {
        auto config_path = resolve_path(vm_init.at("config").as<fs::path>(), vm_init);
        parse_config_file(config_path, vm, all);
    }

    vm_init.clear();
    po::store(run(po::command_line_parser(argc, argv).options(all)), vm);
    validate(vm);
    po::notify(vm);

    return vm;
}

class InvalidWorkingDirectory : public UserError
{
    std::string do_where() const override
    {
        return "get_working_directory";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The working directory you specified ";
        ss << path_;
        ss << " does not exist";
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "enter a valid working directory";
    }
    
    fs::path path_;
public:
    InvalidWorkingDirectory(fs::path p) : path_ {std::move(p)} {}
};

fs::path get_working_directory(const OptionMap& options)
{
    if (options.count("working-directory") == 1) {
        auto result = expand_user_path(options.at("working-directory").as<fs::path>());
        if (!fs::exists(result) && !fs::is_directory(result)) {
            throw InvalidWorkingDirectory {result};
        }
        return result;
    }
    return fs::current_path();
}

fs::path resolve_path(const fs::path& path, const OptionMap& options)
{
    return ::octopus::resolve_path(path, get_working_directory(options));
}

namespace {

std::string prepend_dashes(std::string option)
{
    option.insert(0, "--");
    return option;
}

std::string implode(std::vector<std::string> options)
{
    std::transform(std::cbegin(options), std::cend(options), std::begin(options), prepend_dashes);
    static const std::string delim {" | "};
    return utils::join(options, delim);
}

}

class CommandLineError : public UserError
{
public:
    CommandLineError() = default;
    
    CommandLineError(std::string&& why) : why_ {std::move(why)} {}
    
protected:
    std::string why_;

private:
    virtual std::string do_where() const override
    {
        return "parse_options";
    }
    
    virtual std::string do_why() const override
    {
        return why_;
    }
    
    virtual std::string do_help() const override
    {
        return "use the --help command to view required and allowable options";
    }
};

class BadConfigFile : public  CommandLineError
{
public:
    BadConfigFile(fs::path p)
    {
        std::ostringstream ss {};
        ss << "The config file path (" << p << ") given in the option '--config' does not exist";
        why_ = ss.str();
    }
};

class UnknownCommandLineOption : public CommandLineError
{
public:
    UnknownCommandLineOption(std::string option)
    : CommandLineError { "The option you specified '--" + option + "' is not recognised"}
    {}
};

class MissingRequiredCommandLineArguement : public CommandLineError
{
public:
    MissingRequiredCommandLineArguement(std::string option)
    : CommandLineError {"The command line option '--" + option + "' is required but is missing"}
    {}
    
    MissingRequiredCommandLineArguement(std::vector<std::string> options, bool strict = false)
    {
        std::ostringstream ss {};
        if (strict) {
            ss << "One ";
        } else {
            ss << "At least one ";
        }
        ss << "of the command line options '" + implode(options) + "' is required but none are present";
        why_ = ss.str();
    }
};

class InvalidCommandLineOptionValue : public CommandLineError
{
public:
    template <typename T>
    InvalidCommandLineOptionValue(std::string option, T value, std::string reason)
    : CommandLineError {
    "The arguement '" + std::to_string(value) + "' given to option '--" + option
    + "' was rejected as it " + reason
    } {}
};

class ConflictingCommandLineOptions : public CommandLineError
{
public:
    ConflictingCommandLineOptions(std::vector<std::string> conflicts)
    {
        std::ostringstream ss {};
        ss << "the options";
        for (const auto& option : conflicts) {
            ss << " " << option;
        }
        ss << " are mutually exclusive";
        why_ = ss.str();
    }
};

class MissingDependentCommandLineOption : public CommandLineError
{
public:
    MissingDependentCommandLineOption(std::string given, std::string dependent)
    {
        std::ostringstream ss {};
        ss << "The option " << given << " requires option " << dependent;
        why_ = ss.str();
    }
};

void parse_config_file(const fs::path& config_file, OptionMap& vm, const po::options_description& options)
{
    if (!fs::exists(config_file)) {
        throw BadConfigFile {config_file};
    }
    std::ifstream config {config_file.string()};
    if (config) {
        try {
            po::store(po::parse_config_file(config, options), vm);
        } catch (const po::invalid_config_file_syntax& e) {
            throw CommandLineError {e.what()};
        } catch (const po::unknown_option& e) {
            throw UnknownCommandLineOption {po::strip_prefixes(e.get_option_name())};
        } catch (const po::invalid_option_value& e) {
            throw CommandLineError {e.what()};
        } catch (const po::invalid_bool_value& e) {
            throw CommandLineError {e.what()};
        } catch (const po::ambiguous_option& e) {
            throw CommandLineError {e.what()};
        } catch (const po::reading_file& e) {
            throw CommandLineError {e.what()};
        }
    }
}

void check_positive(const std::string& option, const OptionMap& vm)
{
    if (vm.count(option) == 1) {
        const auto value = vm.at(option).as<int>();
        if (value < 0) {
            throw InvalidCommandLineOptionValue {option, value, "must be positive" };
        }
    }
}

void check_strictly_positive(const std::string& option, const OptionMap& vm)
{
    if (vm.count(option) == 1) {
        const auto value = vm.at(option).as<int>();
        if (value < 1) {
            throw InvalidCommandLineOptionValue {option, value, "must be greater than zero" };
        }
    }
}

void check_probability(const std::string& option, const OptionMap& vm)
{
    if (vm.count(option) == 1) {
        const auto value = vm.at(option).as<float>();
        if (value < 0 || value > 1) {
            throw InvalidCommandLineOptionValue {option, value, "must be between zero and one" };
        }
    }
}

void conflicting_options(const OptionMap& vm, const std::string& opt1, const std::string& opt2)
{
    if (vm.count(opt1) == 1 && !vm[opt1].defaulted() && vm.count(opt2) == 1 && !vm[opt2].defaulted()) {
        throw ConflictingCommandLineOptions {{opt1, opt2}};
    }
}

void option_dependency(const OptionMap& vm, const std::string& given, const std::string& dependent)
{
    if (vm.count(given) == 1 && !vm[given].defaulted())
        if (vm.count(dependent) == 0 || vm[dependent].defaulted()) {
            throw MissingDependentCommandLineOption {given, dependent};
        }
}

void check_reads_present(const OptionMap& vm)
{
    if (vm.count("reads") == 0 && vm.count("reads-file") == 0) {
        throw MissingRequiredCommandLineArguement {std::vector<std::string> {"reads", "reads-file"}};
    }
}

void check_region_files_consistent(const OptionMap& vm)
{
    if (vm.count("regions-file") == 1 && vm.count("skip-regions-file") == 1) {
        const auto regions_file = vm.at("regions-file").as<fs::path>();
        const auto skip_regions_file = vm.at("skip-regions-file").as<fs::path>();
        if (regions_file == skip_regions_file) {
            throw std::invalid_argument {"options 'regions-file' and 'skip-regions-file' must be unique"};
        }
    }
}

void check_trio_consistent(const OptionMap& vm)
{
    if (vm.at("caller").as<std::string>() == "trio"
        && (vm.count("maternal-sample") == 0 || vm.count("paternal-sample") == 0)) {
        throw std::logic_error {"option 'maternal-sample' and 'paternal-sample' are required"
            " when caller=trio"};
    }
}

void validate_caller(const OptionMap& vm)
{
    if (vm.count("caller") == 1) {
        const auto caller = vm.at("caller").as<std::string>();
        static const std::array<std::string, 4> validCallers {
            "individual", "population", "cancer", "trio"
        };
        if (std::find(std::cbegin(validCallers), std::cend(validCallers), caller) == std::cend(validCallers)) {
            throw po::validation_error {po::validation_error::kind_t::invalid_option_value, caller,
                "caller"};
        }
    }
}

po::parsed_options run(po::command_line_parser& parser)
{
    try {
        return parser.run();
    } catch (const po::required_option& e) {
        throw MissingRequiredCommandLineArguement {po::strip_prefixes(e.get_option_name())};
    } catch (const po::unknown_option& e) {
        throw UnknownCommandLineOption {po::strip_prefixes(e.get_option_name())};
    } catch (const po::invalid_option_value& e) {
        throw CommandLineError {e.what()};
    } catch (const po::invalid_bool_value& e) {
        throw CommandLineError {e.what()};
    } catch (const po::ambiguous_option& e) {
        throw CommandLineError {e.what()};
    } catch (const po::reading_file& e) {
        throw CommandLineError {e.what()};
    } catch (const po::invalid_command_line_syntax& e) {
        throw CommandLineError {e.what()};
    } catch (const po::error& e) {
        throw CommandLineError {e.what()};
    }
}

void validate(const OptionMap& vm)
{
    const std::vector<std::string> positive_int_options {
        "threads", "mask-low-quality-tails", "soft-clip-mask-threshold", "mask-soft-clipped-boundary-bases",
        "min-mapping-quality", "good-base-quality", "min-good-bases", "min-read-length",
        "max-read-length", "min-base-quality", "min-supporting-reads", "max-variant-size",
        "num-fallback-kmers", "max-assemble-region-overlap", "assembler-mask-base-quality",
        "min-kmer-prune", "max-bubbles", "max-holdout-depth"
    };
    const std::vector<std::string> strictly_positive_int_options {
        "max-open-read-files", "downsample-above", "downsample-target",
        "max-region-to-assemble", "fallback-kmer-gap", "organism-ploidy",
        "max-haplotypes", "haplotype-holdout-threshold", "haplotype-overflow",
        "max-joint-genotypes"
    };
    const std::vector<std::string> probability_options {
        "snp-heterozygosity", "snp-heterozygosity-stdev", "indel-heterozygosity",
        "somatic-mutation-rate", "min-expected-somatic-frequency", "min-credible-somatic-frequency", "credible-mass",
        "snv-denovo-mutation-rate", "indel-denovo-mutation-rate"
    };
    conflicting_options(vm, "maternal-sample", "normal-sample");
    conflicting_options(vm, "paternal-sample", "normal-sample");
    for (const auto& option : positive_int_options) {
        check_positive(option, vm);
    }
    for (const auto& option : strictly_positive_int_options) {
        check_strictly_positive(option, vm);
    }
    for (const auto& option : probability_options) {
        check_probability(option, vm);
    }
    check_reads_present(vm);
    check_region_files_consistent(vm);
    check_trio_consistent(vm);
    validate_caller(vm);
}

std::istream& operator>>(std::istream& in, ContigPloidy& result)
{
    static const std::regex re {"(?:([^:]*):)?([^=]+)=(\\d+)"};
    
    std::string token;
    in >> token;
    std::smatch match;
    
    if (std::regex_match(token, match, re) && match.size() == 4) {
        if (match.length(1) > 0) {
            result.sample = match.str(1);
        }
        result.contig = match.str(2);
        result.ploidy = boost::lexical_cast<decltype(result.ploidy)>(match.str(3));
    } else {
        using Error = po::validation_error;
        throw Error {Error::kind_t::invalid_option_value, token, "contig-ploidies"};
    }
    
    return in;
}

std::ostream& operator<<(std::ostream& out, const ContigPloidy& cp)
{
    if (cp.sample) out << *cp.sample << ':';
    out << cp.contig << "=" << cp.ploidy;
    return out;
}

std::istream& operator>>(std::istream& in, RefCallType& result)
{
    std::string token;
    in >> token;
    if (token == "positional")
        result = RefCallType::positional;
    else if (token == "blocked")
        result = RefCallType::blocked;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "refcalls"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const RefCallType& type)
{
    switch (type) {
        case RefCallType::positional:
            out << "positional";
            break;
        case RefCallType::blocked:
            out << "blocked";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ContigOutputOrder& result)
{
    std::string token;
    in >> token;
    if (token == "lexicographicalAscending")
        result = ContigOutputOrder::lexicographicalAscending;
    else if (token == "lexicographicalDescending")
        result = ContigOutputOrder::lexicographicalDescending;
    else if (token == "contigSizeAscending")
        result = ContigOutputOrder::contigSizeAscending;
    else if (token == "contigSizeDescending")
        result = ContigOutputOrder::contigSizeDescending;
    else if (token == "asInReference")
        result = ContigOutputOrder::asInReferenceIndex;
    else if (token == "asInReferenceReversed")
        result = ContigOutputOrder::asInReferenceIndexReversed;
    else if (token == "unspecified")
        result = ContigOutputOrder::unspecified;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "contig-output-order"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ContigOutputOrder& order)
{
    switch (order) {
        case ContigOutputOrder::lexicographicalAscending:
            out << "lexicographicalAscending";
            break;
        case ContigOutputOrder::lexicographicalDescending:
            out << "lexicographicalDescending";
            break;
        case ContigOutputOrder::contigSizeAscending:
            out << "contigSizeAscending";
            break;
        case ContigOutputOrder::contigSizeDescending:
            out << "contigSizeDescending";
            break;
        case ContigOutputOrder::asInReferenceIndex:
            out << "asInReferenceIndex";
            break;
        case ContigOutputOrder::asInReferenceIndexReversed:
            out << "asInReferenceIndexReversed";
            break;
        case ContigOutputOrder::unspecified:
            out << "unspecified";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ExtensionLevel& result)
{
    std::string token;
    in >> token;
    if (token == "conservative")
        result = ExtensionLevel::conservative;
    else if (token == "normal")
        result = ExtensionLevel::normal;
    else if (token == "optimistic")
        result = ExtensionLevel::optimistic;
    else if (token == "aggressive")
        result = ExtensionLevel::aggressive;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "extension-level"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ExtensionLevel& level)
{
    switch (level) {
    case ExtensionLevel::conservative:
        out << "conservative";
        break;
    case ExtensionLevel::normal:
        out << "normal";
        break;
    case ExtensionLevel::optimistic:
        out << "optimistic";
        break;
    case ExtensionLevel::aggressive:
        out << "aggressive";
        break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, PhasingLevel& result)
{
    std::string token;
    in >> token;
    if (token == "minimal")
        result = PhasingLevel::minimal;
    else if (token == "conservative")
        result = PhasingLevel::conservative;
    else if (token == "moderate")
        result = PhasingLevel::moderate;
    else if (token == "normal")
        result = PhasingLevel::normal;
    else if (token == "aggressive")
        result = PhasingLevel::aggressive;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "phasing-level"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const PhasingLevel& level)
{
    switch (level) {
        case PhasingLevel::minimal:
            out << "minimal";
            break;
        case PhasingLevel::conservative:
            out << "conservative";
            break;
        case PhasingLevel::moderate:
            out << "moderate";
            break;
        case PhasingLevel::normal:
            out << "normal";
            break;
        case PhasingLevel::aggressive:
            out << "aggressive";
            break;
    }
    return out;
}

} // namespace options
} // namespace octopus
