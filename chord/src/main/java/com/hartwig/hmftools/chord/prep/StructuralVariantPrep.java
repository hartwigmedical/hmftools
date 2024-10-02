package com.hartwig.hmftools.chord.prep;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;

import java.nio.file.NoSuchFileException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.variant.StructuralVariant;
import com.hartwig.hmftools.chord.variant.VcfFile;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.SvVcfTags;

import htsjdk.variant.variantcontext.VariantContext;

public class StructuralVariantPrep
{
    private final ChordConfig mConfig;

    //private static final List<StructuralVariantType> SV_TYPES = List.of(DEL, DUP, INV, BND);
    private static final int[] SV_LENGTH_INTERVALS = { 0, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, Integer.MAX_VALUE};
    private static final DecimalFormat SV_LENGTH_FORMAT = new DecimalFormat("0E00");

    private final Map<SvTypeLengthBin, Integer> mTypeLengthCounts = initializeSizeBins();

    private static final SvTypeLengthBin TRANSLOCATION_BIN = new SvTypeLengthBin(BND, Integer.MIN_VALUE, Integer.MIN_VALUE);
    private static final String TRANSLOCATION = "TRA";

    public StructuralVariantPrep(ChordConfig config)
    {
        mConfig = config;
    }

    public List<StructuralVariant> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.purpleSvVcfFile(sampleId), true);
        List<VariantContext> variantContexts = vcfFile.loadVariants();

        // Only keep the first mate of breakend pairs
        Set<String> idSet = new HashSet<>();
        List<StructuralVariant> variants = new ArrayList<>();
        for(VariantContext variantContext : variantContexts)
        {
            String id = variantContext.getID();
            String mateId = variantContext.getAttributeAsString(SvVcfTags.MATE_ID, id); // For SGLs, use id as mate id

            if(!idSet.contains(id))
            {
                variants.add(new StructuralVariant(variantContext));
                CHORD_LOGGER.trace("  Added SV with id({}) mateId({}): {} ", id, mateId, variantContext);
            }
            else
            {
                CHORD_LOGGER.trace("Skipped SV with id({}) mateId({}): {} ", id, mateId, variantContext);
            }

            idSet.add(id);
            idSet.add(mateId);
        }

        CHORD_LOGGER.debug("Loaded {} structural variants", variants.size());

        return variants;
    }

    public List<MutTypeCount> extractSampleData(String sampleId)
    {
        try
        {
            List<StructuralVariant> variants = loadVariants(sampleId);

            for(StructuralVariant variant : variants)
            {
                StructuralVariantType type = variant.Type;

                if(type == DEL || type == DUP || type == INV || type == BND)
                {
                    SvTypeLengthBin bin = getSvTypeLengthBin(type, variant.RefPosition, variant.AltPosition);

                    CHORD_LOGGER.trace("{} bin({})", variant, bin);

                    mTypeLengthCounts.compute(bin, (k,v) -> v + 1);
                }
            }

            List<MutTypeCount> counts = new ArrayList<>();

            for (SvTypeLengthBin bin: mTypeLengthCounts.keySet()) {
                String binString = bin.toString();
                int count = mTypeLengthCounts.get(bin);

                counts.add(new MutTypeCount(binString, count));

                CHORD_LOGGER.trace("{} {}", binString, String.valueOf(count));
            }

            return counts;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("sample({}) failed to extract SV contexts:", sampleId);
            e.printStackTrace();
            System.exit(1);
            return null;
        }
    }

    private static class SvTypeLengthBin
    {
        StructuralVariantType Type;
        int LowerInterval;
        int UpperInterval;

        SvTypeLengthBin(StructuralVariantType type, int lowerInterval, int upperInterval)
        {
            Type = type;
            LowerInterval = lowerInterval;
            UpperInterval = upperInterval;
        }

        @Override
        public String toString()
        {
            if(Type == BND)
                return TRANSLOCATION;

            String lowerIntervalString = SV_LENGTH_FORMAT.format(LowerInterval);
            String upperIntervalString = (UpperInterval != Integer.MAX_VALUE) ? SV_LENGTH_FORMAT.format(UpperInterval) : "Inf";

            return String.join("_", Type.toString(), lowerIntervalString, upperIntervalString, "bp");
        }
    }

    public static Map<SvTypeLengthBin, Integer> initializeSizeBins()
    {
        Map<SvTypeLengthBin, Integer> bins = new LinkedHashMap<>();

        for(StructuralVariantType type : List.of(DEL, DUP, INV))
        {
            for(int i = 0; i < SV_LENGTH_INTERVALS.length-1; i++)
            {
                int lowerInterval = SV_LENGTH_INTERVALS[i];
                int upperInterval = SV_LENGTH_INTERVALS[i + 1];

                SvTypeLengthBin bin = new SvTypeLengthBin(type, lowerInterval, upperInterval);

                bins.put(bin, 0);
            }
        }

        bins.put(TRANSLOCATION_BIN, 0);

        return bins;
    }

    public SvTypeLengthBin getSvTypeLengthBin(StructuralVariantType svType, int startPos, int endPos)
    {
        if(svType == BND)
            return TRANSLOCATION_BIN;

        int svLength = Math.abs(startPos - endPos);

        for(SvTypeLengthBin bin : mTypeLengthCounts.keySet())
        {
            if(svType == bin.Type && svLength > bin.LowerInterval && svLength <= bin.UpperInterval)
                return bin;
        }

        throw new IllegalArgumentException(String.format("Failed to get SV type/length bin for svType(%s) svLength(%s)", svType, svLength));
    }
}
