package com.hartwig.hmftools.esvee.utils.vcfcompare;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SvVcfComparer
{
    private final String mSampleId;
    private final String mReferenceId;

    private final String mOldVcf;
    private final String mNewVcf;

    private final String mOldUnfilteredVcf;
    private final String mNewUnfilteredVcf;

    private final String mOutputDir;
    private final String mOutputId;

    private final Map<String,List<VariantBreakend>> mOldChrBreakendMap;
    private final Map<String,List<VariantBreakend>> mNewChrBreakendMap;

    private final Map<String,List<VariantBreakend>> mOldChrBreakendMapUnfiltered;
    private final Map<String,List<VariantBreakend>> mNewChrBreakendMapUnfiltered;

    private final boolean mShowNonPass;

    private final BufferedWriter mWriter;

    private RefGenomeVersion mRefGenomeVersion;

    private static final String OLD_VCF = "old_vcf";
    private static final String NEW_VCF = "new_vcf";

    private static final String OLD_UNFILTERED_VCF = "old_unfiltered_vcf";
    private static final String NEW_UNFILTERED_VCF = "new_unfiltered_vcf";

    private static final String SHOW_NON_PASS = "show_non_pass";

    private static final int DEFAULT_MAX_DIFF = 20;
    private static final double DEFAULT_MAX_DIFF_PERC = 0.2;
    private static final int DEFAULT_APPROX_MATCH_UPPER_LOWER_BOUNDS = 10;

    public SvVcfComparer(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mReferenceId = configBuilder.getValue(REFERENCE, "");

        mOldVcf = configBuilder.getValue(OLD_VCF);
        mNewVcf = configBuilder.getValue(NEW_VCF);

        mOldUnfilteredVcf = configBuilder.getValue(OLD_UNFILTERED_VCF);
        mNewUnfilteredVcf = configBuilder.getValue(NEW_UNFILTERED_VCF);

        mShowNonPass = configBuilder.hasFlag(SHOW_NON_PASS);

        mOutputDir = FileWriterUtils.parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mOldChrBreakendMap = new HashMap<>();
        mNewChrBreakendMap = new HashMap<>();

        mOldChrBreakendMapUnfiltered = new HashMap<>();
        mNewChrBreakendMapUnfiltered = new HashMap<>();

        mRefGenomeVersion = null;

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(mOldVcf == null || mNewVcf == null)
        {
            SV_LOGGER.error("Missing VCFs");
            return;
        }

        loadVariants(mOldVcf, mOldChrBreakendMap);
        loadVariants(mNewVcf, mNewChrBreakendMap);

        matchAndWriteVariants(mOldChrBreakendMap, mNewChrBreakendMap);

        if(mOldUnfilteredVcf != null)
        {
            loadVariants(mOldUnfilteredVcf, mOldChrBreakendMapUnfiltered);
            matchAndWriteVariants(mNewChrBreakendMap, mOldChrBreakendMapUnfiltered);
        }

        if(mNewUnfilteredVcf != null)
        {
            loadVariants(mNewUnfilteredVcf, mNewChrBreakendMapUnfiltered);
            matchAndWriteVariants(mOldChrBreakendMap, mNewChrBreakendMapUnfiltered);
        }

        writeUnmatchedVariants(mOldChrBreakendMap, true);
        writeUnmatchedVariants(mNewChrBreakendMap, false);

        FileWriterUtils.closeBufferedWriter(mWriter);

        SV_LOGGER.info("Esvee compare VCFs complete");
    }

    private void loadVariants(final String vcfFile, final Map<String,List<VariantBreakend>> chrBreakendMap)
    {
        SV_LOGGER.info("Loading vcfFile({})", vcfFile);

        VcfFileReader reader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = reader.vcfHeader();
        GenotypeIds genotypeIds = GenotypeIds.fromVcfHeader(vcfHeader, mReferenceId, mSampleId);

        if(genotypeIds == null)
        {
            System.exit(1);
        }

        SV_LOGGER.info("Genotype info: ref({}: {}) tumor({}: {})",
                genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);

        String currentChr = "";
        List<VariantBreakend> breakends = null;

        SvCaller svCaller = SvCaller.fromVcfPath(vcfFile);
        VcfType sourceVcfType = VcfType.fromVcfPath(vcfFile);

        for(VariantContext variantContext : reader.iterator())
        {
            String chromosome = variantContext.getContig();

            if(mRefGenomeVersion == null)
                mRefGenomeVersion = chromosome.startsWith(CHR_PREFIX) ? V38 : V37;

            if(!currentChr.equals(chromosome))
            {
                currentChr = chromosome;
                breakends = new ArrayList<>();
                chrBreakendMap.put(chromosome, breakends);
            }

            breakends.add(new VariantBreakend(variantContext, svCaller, sourceVcfType));
        }

        SV_LOGGER.info("Loaded {} SVs from {})",
                chrBreakendMap.values().stream().mapToInt(x -> x.size()).sum(), vcfFile);
    }

    private void matchAndWriteVariants(
            Map<String,List<VariantBreakend>> oldChrBreakendMap,
            Map<String,List<VariantBreakend>> newChrBreakendMap,
            MatchType matchType
    )
    {
        int matchedCount = 0;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

            List<VariantBreakend> oldBreakends = oldChrBreakendMap.get(chrStr);
            List<VariantBreakend> newBreakends = newChrBreakendMap.get(chrStr);

            if(oldBreakends == null || newBreakends == null)
                continue;

            for(int oldBreakendIndex = 0; oldBreakendIndex < oldBreakends.size(); oldBreakendIndex++)
            {
                VariantBreakend oldBreakend = oldBreakends.get(oldBreakendIndex);

                for(int newBreakendIndex = 0; newBreakendIndex < newBreakends.size(); newBreakendIndex++)
                {
                    VariantBreakend newBreakend = newBreakends.get(newBreakendIndex);

                    boolean hasMatch = false;

                    switch(matchType)
                    {
                        case EXACT_MATCH:
                            hasMatch = oldBreakend.exactMatch(newBreakend);
                            break;
                        case COORDS_ONLY:
                            hasMatch = oldBreakend.coordsOnlyMatch(newBreakend);
                            break;
                        case APPROX_MATCH:
                            hasMatch = oldBreakend.approxMatch(newBreakend, DEFAULT_APPROX_MATCH_UPPER_LOWER_BOUNDS);
                            break;
                    }

                    if(hasMatch)
                    {
                        if(mShowNonPass || oldBreakend.isPassVariant() || newBreakend.isPassVariant())
                            writeBreakend(oldBreakend, newBreakend, matchType);

                        // Only allow one match
                        oldBreakends.remove(oldBreakendIndex);
                        newBreakends.remove(newBreakendIndex);

                        // Decrement index to account for removal of breakends
                        oldBreakendIndex--;
                        newBreakendIndex--;

                        matchedCount++;

                        break;
                    }
                }
            }
        }

        if(matchedCount > 0)
        {
            SV_LOGGER.debug("  Found {} variants with match type: {}", matchedCount, matchType);
        }
    }

    private void  matchAndWriteVariants(
            Map<String,List<VariantBreakend>> oldChrBreakendMap,
            Map<String,List<VariantBreakend>> newChrBreakendMap
    )
    {
        matchAndWriteVariants(oldChrBreakendMap, newChrBreakendMap, MatchType.EXACT_MATCH);
        matchAndWriteVariants(oldChrBreakendMap, newChrBreakendMap, MatchType.COORDS_ONLY);
        matchAndWriteVariants(oldChrBreakendMap, newChrBreakendMap, MatchType.APPROX_MATCH);
    }

    private void writeUnmatchedVariants(Map<String,List<VariantBreakend>> chrBreakendMap, boolean isOld)
    {
        int unmatchedVariantsCount = 0;

        for(List<VariantBreakend> breakends : chrBreakendMap.values())
        {
            if(breakends == null)
                continue;

            for(VariantBreakend breakend : breakends)
            {
                if(mShowNonPass || breakend.isPassVariant())
                {
                    if(isOld)
                        writeBreakend(breakend, null, MatchType.NO_MATCH);
                    else
                        writeBreakend(null, breakend, MatchType.NO_MATCH);
                }

                unmatchedVariantsCount++;
            }
        }

        if(unmatchedVariantsCount > 0)
        {
            SV_LOGGER.debug("Writing {} unmatched variants", unmatchedVariantsCount);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + ".sv_compare";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += TSV_EXTENSION;

            SV_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = FileWriterUtils.createBufferedWriter(fileName, false);

            String header = String.join(
                    TSV_DELIM, "SampleId",
                    "OldId",       "NewId",
                    "MatchType",
                    "Diffs",
                    "OldSvCoords", "NewSvCoords",
                    "OldCoords",   "NewCoords",
                    "OldCipos",    "NewCipos",
                    "OldIhompos",  "NewIhompos",
                    "OldHomSeq",   "NewHomSeq",
                    "OldInsSeq",   "NewInsSeq",
                    "OldSvType",   "NewSvType",
                    "OldFilter",   "NewFilter",
                    "OldVcfType",  "NewVcfType",
                    "OldQual",     "NewQual",
                    "OldVF",       "NewVF"
            );

            writer.write(header);
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeBreakend(VariantBreakend oldBreakend, VariantBreakend newBreakend, MatchType matchType)
    {
        String oldId       = "";
        String oldSvCoords = "";
        String oldCoords   = "";
        String oldCipos    = "";
        String oldIhompos  = "";
        String oldHomSeq   = "";
        String oldInsSeq   = "";
        String oldSvtype   = "";
        String oldFilter   = "";
        String oldVcfType  = "";
        String oldQual     = "";
        String oldVF       = "";
        if(oldBreakend != null)
        {
            oldId = oldBreakend.Context.getID();
            oldSvCoords = oldBreakend.svCoordStr();
            oldCoords = oldBreakend.coordStr();
            oldCipos = Arrays.toString(oldBreakend.Cipos);
            oldIhompos = Arrays.toString(oldBreakend.Ihompos);
            oldHomSeq = oldBreakend.Homseq;
            oldInsSeq = oldBreakend.AltCoords.InsertSequence;
            oldSvtype = oldBreakend.SvType;
            oldFilter = oldBreakend.filtersStr();
            oldVcfType = oldBreakend.SourceVcfType.toString();
            oldQual = oldBreakend.qualStr();
            oldVF = oldBreakend.getExtendedAttributeAsString(mSampleId, TOTAL_FRAGS);
        }

        String newId       = "";
        String newSvCoords = "";
        String newCoords   = "";
        String newCipos    = "";
        String newIhompos  = "";
        String newHomSeq   = "";
        String newInsSeq   = "";
        String newSvtype   = "";
        String newFilter   = "";
        String newVcfType  = "";
        String newQual     = "";
        String newVF       = "";
        if(newBreakend != null)
        {
            newId = newBreakend.Context.getID();
            newSvCoords = newBreakend.svCoordStr();
            newCoords = newBreakend.coordStr();
            newCipos = Arrays.toString(newBreakend.Cipos);
            newIhompos = Arrays.toString(newBreakend.Ihompos);
            newHomSeq = newBreakend.Homseq;
            newInsSeq = newBreakend.AltCoords.InsertSequence;
            newSvtype = newBreakend.SvType;
            newFilter = newBreakend.filtersStr();
            newVcfType = newBreakend.SourceVcfType.toString();
            newQual = newBreakend.qualStr();
            newVF = newBreakend.getExtendedAttributeAsString(mSampleId, TOTAL_FRAGS);
        }

        String diffs = "";
        if(oldBreakend != null && newBreakend != null)
        {
            diffs = String.join(",", compareBreakendAttributes(oldBreakend, newBreakend, matchType));
        }

        try
        {
            String line = String.join(
                    TSV_DELIM,
                    mSampleId,
                    oldId,      newId,
                    matchType.toString(),
                    diffs,
                    oldSvCoords, newSvCoords,
                    oldCoords,   newCoords,
                    oldCipos,    newCipos,
                    oldIhompos,  newIhompos,
                    oldHomSeq,   newHomSeq,
                    oldInsSeq,   newInsSeq,
                    oldSvtype,   newSvtype,
                    oldFilter,   newFilter,
                    oldVcfType,  newVcfType,
                    oldQual,     newQual,
                    oldVF,       newVF
            );

            mWriter.write(line);
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private static boolean hasDiffWithinTolerance(double value1, double value2, double maxDiff, double maxDiffPerc)
    {
        double diff = abs(value1 - value2);
        double diffPerc = diff / max(value1, value2);
        return diff > maxDiff && diffPerc > maxDiffPerc;
    }

    private static final String DIFF_PASS = "PASS_FILTER";
    private static final String DIFF_COORDS = "COORDS";
    private static final String DIFF_INSSEQ = "INSSEQ";

    private List<String> compareBreakendAttributes(VariantBreakend oldBreakend, VariantBreakend newBreakend, MatchType matchType)
    {
        List<String> diffSet = new ArrayList<>();

        if((oldBreakend.isPassVariant() & !newBreakend.isPassVariant()) || (!oldBreakend.isPassVariant() & newBreakend.isPassVariant()))
        {
            diffSet.add(DIFF_PASS);
        }

        if((matchType == MatchType.APPROX_MATCH ||  matchType == MatchType.COORDS_ONLY))
        {
            if(!oldBreakend.coordStr().equals(newBreakend.coordStr()))
                diffSet.add(DIFF_COORDS);

            if(!Arrays.equals(oldBreakend.Cipos, newBreakend.Cipos))
                diffSet.add(CIPOS);

            if(!Arrays.equals(oldBreakend.Ihompos, newBreakend.Ihompos))
                diffSet.add(IHOMPOS);

            if(!oldBreakend.Homseq.equals(newBreakend.Homseq))
                diffSet.add(HOMSEQ);

            if(!oldBreakend.InsertSequence.equals(newBreakend.InsertSequence))
                diffSet.add(DIFF_INSSEQ);

            if(!oldBreakend.SvType.equals(newBreakend.SvType))
                diffSet.add(SVTYPE);
        }

        if(hasDiffWithinTolerance(
                oldBreakend.getExtendedAttributeAsDouble(mSampleId, TOTAL_FRAGS),
                newBreakend.getExtendedAttributeAsDouble(mSampleId, TOTAL_FRAGS),
                DEFAULT_MAX_DIFF,
                DEFAULT_MAX_DIFF_PERC
        ))
        {
            diffSet.add(TOTAL_FRAGS);
        }

        return diffSet;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);

        configBuilder.addPath(OLD_VCF, true, "Path to the old VCF file");
        configBuilder.addPath(NEW_VCF, true, "Path to the new VCF file");

        configBuilder.addPath(OLD_UNFILTERED_VCF, false, "Path to the old unfiltered VCF file");
        configBuilder.addPath(NEW_UNFILTERED_VCF, false, "Path to the new unfiltered VCF file");

        configBuilder.addPath(SHOW_NON_PASS, false, "Show variants not PASSing in both old nor new VCF files");

        FileWriterUtils.addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        ConfigUtils.setLogLevel(configBuilder);

        SvVcfComparer svVcfCompare = new SvVcfComparer(configBuilder);
        svVcfCompare.run();
    }
}