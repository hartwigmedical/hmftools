package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ReferenceMutation extends Variant
{
    private final String mChromosome;
    private final int mPosition;
    private final String mRef;
    private final String mAlt;
    private final String mSource;
    private final String mGene;

    private static final Logger LOGGER = LogManager.getLogger(ReferenceMutation.class);

    public ReferenceMutation(
            final String chromosome, final int position, final String ref, final String alt, final String source, final String gene)
    {
        mChromosome = chromosome;
        mPosition = position;
        mRef = ref;
        mAlt = alt;
        mSource = source;
        mGene = gene;
    }

    @Override
    public CategoryType categoryType()
    {
        return CategoryType.REFERENCE;
    }

    @Override
    public String description()
    {
        if(!mRef.isEmpty() && !mAlt.isEmpty())
        {
            return format("%s:%s %s>%s %s", mChromosome, mPosition, mRef, mAlt, mSource);
        }
        else
        {
            return format("%s:%s %s", mChromosome, mPosition, mSource);
        }
    }

    @Override
    public String gene()
    {
        return mGene;
    }

    @Override
    public double copyNumber()
    {
        return 0;
    }

    @Override
    public double vaf()
    {
        return 0;
    }

    @Override
    public int tumorFragments()
    {
        return 0;
    }

    @Override
    public boolean reported()
    {
        return false;
    }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome)
    {
        if(mAlt.isEmpty())
        {
            int startLength = PROBE_LENGTH / 2;
            int startPos = mPosition - startLength;

            setSequence(refGenome.getBaseString(mChromosome, startPos, startPos + PROBE_LENGTH - 1));
            return;
        }

        int altLength = mAlt.length();
        int refLength = mRef.length();
        int startLength = PROBE_LENGTH / 2 - altLength / 2;
        int startPos = mPosition - startLength;

        String basesStart = refGenome.getBaseString(mChromosome, startPos, mPosition - 1);
        int endBaseLength = PROBE_LENGTH - basesStart.length() - altLength;

        int postPosition = mPosition + refLength;
        String basesEnd = refGenome.getBaseString(mChromosome, postPosition, postPosition + endBaseLength - 1);

        setSequence(basesStart + mAlt + basesEnd);
    }

    @Override
    public boolean checkFilters()
    {
        return false;
    }

    @Override
    public boolean checkAndRegisterLocation(final ProximateLocations registeredLocations)
    {
        registeredLocations.addRegisteredLocation(mChromosome, mPosition);
        return true;
    }

    public String toString()
    {
        return description();
    }

    public static List<Variant> loadKnownMutations(final String filename)
    {
        List<Variant> variants = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);

            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int refIndex = fieldsIndexMap.get(FLD_REF);
            int altIndex = fieldsIndexMap.get(FLD_ALT);
            int sourceIndex = fieldsIndexMap.get("Source");
            // int geneIndex = fieldsIndexMap.get("Gene");

            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                variants.add(new ReferenceMutation(
                        values[chrIndex], Integer.parseInt(values[posIndex]), values[refIndex], values[altIndex], values[sourceIndex],
                        ""));
            }

            LOGGER.info("loaded {} reference variants from file({})", variants.size(), filename);
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load reference variants file({}): {}", filename, e.toString());
            return null;
        }

        return variants;
    }
}
