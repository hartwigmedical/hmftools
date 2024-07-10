package com.hartwig.hmftools.esvee.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getByteValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

public class TruthsetAnnotation
{
    private final Map<String,List<TruthsetBreakend>> mChrBreakendMap;

    private static final String TRUTHSET_BREAKENDS = "truthset_breakends";

    public TruthsetAnnotation(final String filename)
    {
        mChrBreakendMap = Maps.newHashMap();

        if(filename != null)
            loadFile(filename);
    }

    public boolean enabled() { return !mChrBreakendMap.isEmpty(); }

    public String findTruthsetAnnotation(final Breakend breakend)
    {
        return findAnnotation(breakend.Chromosome, breakend.Position, breakend.Orient);
    }

    private String findAnnotation(final String chromosome, final int position, final Orientation orientation)
    {
        List<TruthsetBreakend> breakends = mChrBreakendMap.get(chromosome);

        if(breakends == null)
            return NO_TRUTHSET_MATCH;

        TruthsetBreakend match = breakends.stream().filter(x -> x.matches(position, orientation)).findFirst().orElse(null);

        return match != null ? match.asTsv() : NO_TRUTHSET_MATCH;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(TRUTHSET_BREAKENDS, false, "Truthset breakend file for annotation");
    }

    public static String filename(final ConfigBuilder configBuilder) { return configBuilder.getValue(TRUTHSET_BREAKENDS); }

    private void loadFile(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            String header = lines.get(0);
            lines.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);

            for(String line : lines)
            {
                String[] values = line.split(CSV_DELIM, -1);
                String chromosome = values[chrIndex];

                List<TruthsetBreakend> breakends = mChrBreakendMap.get(chromosome);

                if(breakends == null)
                {
                    breakends = Lists.newArrayList();
                    mChrBreakendMap.put(chromosome, breakends);
                }

                String[] ciposString = values[fieldsIndexMap.get("Cipos")].split(ITEM_DELIM,2);

                int[] cipos = new int[] { Integer.parseInt(ciposString[0]), Integer.parseInt(ciposString[1]) };

                String fragsStr = values[fieldsIndexMap.get("VF")];
                int fragments = !fragsStr.isEmpty() && !fragsStr.equals("NA") ? Integer.parseInt(fragsStr) : 0;

                TruthsetBreakend truthsetBreakend = new TruthsetBreakend(
                        chromosome,
                        getIntValue(fieldsIndexMap, "JunctionPosition", values),
                        Orientation.fromByte(getByteValue(fieldsIndexMap, "JunctionOrientation", values)),
                        fragments,
                        cipos,
                        StructuralVariantType.valueOf(values[fieldsIndexMap.get("SvType")]));

                breakends.add(truthsetBreakend);
            }

            SV_LOGGER.info("loaded {} truthset breakends from file: {}",
                    mChrBreakendMap.values().stream().mapToInt(x -> x.size()).sum(), filename);

        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read truthset breakends file({})", filename, exception.toString());
        }
    }

    public static String tsvHeader()
    {
        return "TruthPosition\tTruthFrags\tTruthSvType";
    }

    public static String NO_TRUTHSET_MATCH = "NA\t0\tNA";

    private class TruthsetBreakend
    {
        public final String Chromosome;
        public final int Position;
        public final Orientation Orient;
        public final int Fragments;
        public final int[] Cipos;
        public final StructuralVariantType Type;

        public TruthsetBreakend(
                final String chromosome, final int position, final Orientation orientation, final int fragments, final int[] cipos,
                final StructuralVariantType type)
        {
            Chromosome = chromosome;
            Position = position;
            Orient = orientation;
            Fragments = fragments;
            Cipos = cipos;
            Type = type;
        }

        public boolean matches(final int position, final Orientation orientation)
        {
            if(Orient != orientation)
                return false;

            return positionWithin(position, Position + Cipos[0], Position + Cipos[1]);
        }

        public String asTsv()
        {
            return format("%d\t%d\t%s", Position, Fragments, Type);
        }
    }
}
