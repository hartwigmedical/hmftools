package com.hartwig.hmftools.common.variant;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;

public class SimpleVariant extends BasePosition
{
    public final String Ref;
    public final String Alt;
    public final VariantType Type;

    private final int mIndelLength; // number of bases inserted or deleted, negative for DELs

    public SimpleVariant(final String chromosome, final int position, final String ref, final String alt)
    {
        super(chromosome, position);
        Ref = ref;
        Alt = alt;

        if(ref.length() == alt.length())
        {
            Type = ref.length() == 1 ? VariantType.SNP : VariantType.MNP;
            mIndelLength = 0;
        }
        else
        {
            Type = VariantType.INDEL;
            mIndelLength = alt().length() - ref().length();
        }
    }

    public String chromosome() { return Chromosome; }
    public int position() { return Position; }

    public int positionEnd()
    {
        // the right-aligned positional end of the variant
        if(isSNV())
            return position();
        else if(isInsert())
            return position() + 1;
        else // deletes and MNVs
            return position() + Ref.length() - 1;
    }

    // convenience
    public String ref() { return Ref; }
    public String alt() { return Alt; }

    public int end() { return position() + ref().length() - 1; }

    public boolean isSNV() { return Type == VariantType.SNP; }
    public boolean isMNV() { return Type == VariantType.MNP; }
    public boolean isIndel() { return Type == VariantType.INDEL; }

    public boolean isDelete() { return ref().length() > alt().length(); }
    public boolean isInsert() { return alt().length() > ref().length(); }

    public int altLength() { return Alt.length(); }
    public int refLength() { return Ref.length(); }
    public int indelLength() { return mIndelLength; }
    public int indelLengthAbs() { return abs(mIndelLength); }

    public boolean matches(final SimpleVariant variant) { return matches(variant.Chromosome, variant.Position, variant.Ref, variant.Alt); }

    public boolean matches(final String chromosome, final int position, final String ref, final String alt)
    {
        return Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt);
    }

    public boolean matches(final VariantContextDecorator somaticVariant)
    {
        return somaticVariant.chromosome().equals(Chromosome) && somaticVariant.position() == Position
                && somaticVariant.ref().equals(Ref) && somaticVariant.alt().equals(Alt);
    }

    public String toString() { return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt); }

    private static final int VAR_ITEM_COUNT = 4;
    private static final String VAR_ITEM_DELIM = ":";

    public static List<SimpleVariant> fromConfig(final String configStr)
    {
        String[] variantStrs = configStr.split(ITEM_DELIM, -1);

        List<SimpleVariant> variants = Lists.newArrayList();

        for(String variantStr : variantStrs)
        {
            String[] items = variantStr.split(VAR_ITEM_DELIM, VAR_ITEM_COUNT);

            if(items.length != VAR_ITEM_COUNT)
                return null;

            variants.add(new SimpleVariant(items[0], Integer.parseInt(items[1]), items[2], items[3]));
        }

        return variants;
    }

    public static List<SimpleVariant> loadSimpleVariants(final String filename) throws Exception
    {
        List<SimpleVariant> sampleVariants = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int refIndex = fieldsIndexMap.get(FLD_REF);
            int altIndex = fieldsIndexMap.get(FLD_ALT);

            for(String line : lines)
            {
                String[] values = line.split(CSV_DELIM, -1);

                String ref = values[refIndex];
                String alt = values[altIndex];

                sampleVariants.add(new SimpleVariant(values[chrIndex], Integer.parseInt(values[posIndex]), ref, alt));
            }
        }
        catch(Exception e)
        {
            e.printStackTrace();
            throw new Exception(format("failed to load variants file(%s)", filename));
        }

        return sampleVariants;
    }
}