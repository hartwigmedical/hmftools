package com.hartwig.hmftools.esvee;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.config.Advanced;
import com.hartwig.hmftools.esvee.config.CommandLine;
import com.hartwig.hmftools.esvee.config.HMFConfig;
import com.hartwig.hmftools.esvee.config.Optional;
import org.immutables.value.Value;
import org.jetbrains.annotations.Nullable;

import java.io.File;

@Value.Immutable
public interface SVAConfig extends HMFConfig {
    @CommandLine(description = "Tumor BAM file location")
    File bamFile();

    @CommandLine(name = "germline_bam", description = "Germline BAM file location")
    @Optional(defaultValue = "null")
    @Nullable
    File germlineBAMFile();

    @CommandLine(name = "ref_genome", description = "Reference genome FASTA file")
    File referenceGenomeFile();

    @CommandLine(name = "ref_genome_index", description = "Reference genome index file")
    File referenceGenomeIndex();

    @CommandLine(name = "alt_ref_genome_index", description = "Alternate reference genome index file")
    @Optional(defaultValue = "null")
    @Nullable
    File altReferenceGenomeIndex();

    @CommandLine(description = "CSV containing candidates")
    File junctionFile();

    @CommandLine(description = "Output VCF file")
    String outputFile();

    @Nullable
    @CommandLine(description = "Output CSV file")
    @Optional(defaultValue = "null")
    @Advanced
    String outputCSVFile();

    @Nullable
    @CommandLine(description = "Output BAM file")
    @Optional(defaultValue = "null")
    String outputBAMFile();

    @Optional(defaultValue = "V37")
    RefGenomeVersion referenceGenomeVersion();

    @CommandLine(description = "Whether to create HTML files describing the processing of each junction")
    @Optional(defaultValue = "false")
    boolean createHTMLSummaries();

    @CommandLine(description = "If creating HTML summary files, the maximum number of files that will be produced")
    @Optional(defaultValue = "10000")
    @Advanced
    int maxHTMLSummaries();

    @CommandLine(description = "The folder in which to place HTML summaries")
    @Optional(defaultValue = "Summaries")
    @Advanced
    String htmlSummariesFolder();

    @CommandLine(description = "Whether created HTML files contain SVGs (no effect without createHTMLSummaries)")
    @Optional(defaultValue = "true")
    boolean createDiagrams();

    @CommandLine(description = "How many threads to use for processing. -1 means \"all cores\"")
    @Optional(defaultValue = "1")
    int threads();

    @CommandLine(description = "When a base is considered in the context of a single read, at what quality level do we start to see this base as low-quality")
    @Optional(defaultValue = "26")
    @Advanced
    int lowBaseQualThreshold();

    @CommandLine(description = "When a base has more than 1 piece of supporting evidence, what is the cumulative total below which the base is low quality anyway")
    @Optional(defaultValue = "39")
    @Advanced
    int lowBaseQualCumulativeThreshold();

    @CommandLine(description = "The average baseq below which we consider the entire read (or section of a read) to be low quality")
    @Optional(defaultValue = "30")
    @Advanced
    int averageQualityThreshold();

    @Optional(defaultValue = "1")
    @Advanced
    int maxMismatchedCountForStrongSupport();

    @Optional(defaultValue = "4")
    @Advanced
    int maxMismatchedCountForWeakSupport();

    @Optional(defaultValue = "5")
    @Advanced
    int maxMismatchedCountForDedupingAssemblies();

    @Optional(defaultValue = "50")
    @Advanced
    int maxDistanceToDedupeAssemblies();

    @Optional(defaultValue = "true")
    @Advanced
    boolean tryExtendingUsingDiscordantReads();

    @Optional(defaultValue = "1000")
    @Advanced
    int discordantPairFragmentLength();

    @Optional(defaultValue = "400")
    @Advanced
    int discordantPairSearchDistance();

    @Optional(defaultValue = "31")
    @Advanced
    int discordantPairMinMapQ();

    @Optional(defaultValue = "20")
    @Advanced
    int assemblyExtensionMinMatchedBases();

    @Optional(defaultValue = "4")
    @Advanced
    int assemblyExtensionMaxRepeatScore();

    @Optional(defaultValue = "1")
    @Advanced
    int assemblyExtensionMaxMismatches();

    @Optional(defaultValue = "1")
    @Advanced
    int maxMismatchesForFolding();

    @CommandLine(description = "If a variant or assembly has less than this many fragments support it, it is dropped.")
    @Optional(defaultValue = "2")
    @Advanced
    int minReadsToSupportAssembly();

    @Optional(defaultValue = "20")
    @Advanced
    int minMapQToStartJunction();

    @CommandLine(description = "When trimming polgG/C, how many Gs/Cs must appear on the edge of a read to be trimmed.")
    @Optional(defaultValue = "4")
    @Advanced
    int normaliserPolyGLength();

    @CommandLine(description = "Indels this size or larger near the edge of a read will be converted to soft-clips.")
    @Optional(defaultValue = "6")
    @Advanced
    int normaliserIndelMinSizeToSoftClip();

    @CommandLine(description = "How close to the edge is \"near\" for an indels.")
    @Optional(defaultValue = "16")
    @Advanced
    int normaliserIndelMaxEdgeDistance();

    @CommandLine(description = "Not a MapQ. Ignore any alignments during extension that come from BWA with a score less than this")
    @Optional(defaultValue = "20")
    @Advanced
    int alignerMinScore();

    @CommandLine(description = "If there is more than one candidate for extending an assembly alignment, ignore any that insert this many more bases than the best")
    @Optional(defaultValue = "16")
    @Advanced
    int alignerExtensionInsertTolerance();

    @CommandLine(description = "We don't attempt to call the aligner if we have less than this many bases")
    @Optional(defaultValue = "20")
    @Advanced
    int alignerMinBases();

    @CommandLine(description = "When extending alignments, if we have a candidate match within this many bases of the existing neighbour, prioritise that alignment")
    @Optional(defaultValue = "2000")
    @Advanced
    int alignerMaxDistanceToConsiderNearby();

    @CommandLine(description = "What is the smallest insertion/deletion size to call")
    @Optional(defaultValue = "32")
    @Advanced
    int callerMinSizeToCall();

    @CommandLine(description = "Whether variants with germline support should be dropped")
    @Optional(defaultValue = "false")
    @Advanced
    boolean dropGermline();

    @CommandLine(description = "The threshold below which a LOW_OVERHANG filter will be applied to the VCF")
    @Optional(defaultValue = "20")
    @Advanced
    int vcfLowOverhangThreshold();

    @CommandLine(description = "The threshold below which a LOW_QUALITY filter will be applied to the VCF")
    @Optional(defaultValue = "40")
    @Advanced
    int vcfLowQualityThreshold();

    @Optional(defaultValue = "false")
    @Advanced
    boolean extendPrimaries();

    @CommandLine(description = "Whether individual operations are timed to prevent slow processing")
    @Optional(defaultValue = "false")
    @Advanced
    boolean timeoutsEnabled();

    @CommandLine(description = "Primary Assembly timeout")
    @Optional(defaultValue = "5000")
    @Advanced
    int primaryAssemblyTimeoutMillis();

    @CommandLine(description = "Assembly Extension timeout")
    @Optional(defaultValue = "5000")
    @Advanced
    int extensionTimeoutMillis();

    @SuppressWarnings("BooleanMethodIsAlwaysInverted")
    @Optional(defaultValue = "false")
    @Advanced
    boolean debug();

    static void addConfig(final ConfigBuilder configBuilder) {
        HMFConfig.addConfig(configBuilder, SVAConfig.class);
    }

    static SVAConfig load(final ConfigBuilder configBuilder) {
        return HMFConfig.load(configBuilder, SVAConfig.class, ImmutableSVAConfig.builder());
    }
}
