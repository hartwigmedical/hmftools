package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;

public class SimpleVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;

    public SimpleVariant(final String chromosome, final int position, final String ref, final String alt, final VariantType type)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
    }

    public boolean matches(final SomaticVariant somaticVariant)
    {
        return somaticVariant.Chromosome.equals(Chromosome) && somaticVariant.Position == Position
                && somaticVariant.Ref.equals(Ref) && somaticVariant.Alt.equals(Alt) && somaticVariant.Type == Type;
    }

    public boolean matches(final VariantContextDecorator somaticVariant)
    {
        return somaticVariant.chromosome().equals(Chromosome) && somaticVariant.position() == Position
                && somaticVariant.ref().equals(Ref) && somaticVariant.alt().equals(Alt);
    }

    public String toString()
    {
        return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt);
    }

    public static List<SimpleVariant> loadSimpleVariants(final String filename)
    {
        List<SimpleVariant> sampleVariants = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int refIndex = fieldsIndexMap.get(FLD_REF);
            int altIndex = fieldsIndexMap.get(FLD_ALT);

            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(CSV_DELIM, -1);

                String ref = values[refIndex];
                String alt = values[altIndex];
                VariantType type = VariantType.type(ref, alt);

                sampleVariants.add(new SimpleVariant(values[chrIndex], Integer.parseInt(values[posIndex]), ref, alt, type));
            }
        }
        catch(Exception e)
        {
            CT_LOGGER.error("failed to load variants file({}): {}", filename, e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        return sampleVariants;
    }
}
