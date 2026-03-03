package com.hartwig.hmftools.common.sigs;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.variant.SmallVariant;

import org.jetbrains.annotations.Nullable;

public class SnvSigUtils
{
    public static final int SNV_TRINUCLEOTIDE_BUCKET_COUNT = 96;

    public static Matrix loadSnvSignatures()
    {
        List<String> signatureNames = Lists.newArrayList();
        return loadSnvSignatures(signatureNames);
    }

    public static Matrix loadSnvSignatures(final List<String> signatureNames)
    {
        String sigDefinitionsFile = "/sigs/snv_cosmic_signatures.csv";

        List<String> sigDefinitionLines = new BufferedReader(new InputStreamReader(
                SnvSigUtils.class.getResourceAsStream(sigDefinitionsFile))).lines().collect(Collectors.toList());

        return loadMatrixDataFile(sigDefinitionLines, signatureNames, Lists.newArrayList(), false);
    }

    public static Map<String,String> loadSnvSignatureEtiologies()
    {
        String sigEtiologiesFile = "/sigs/signatures_etiology.tsv";
        Map<String,String> sigEtiologies = Maps.newHashMap();

        List<String> lines = new BufferedReader(new InputStreamReader(
                SnvSigUtils.class.getResourceAsStream(sigEtiologiesFile))).lines().collect(Collectors.toList());

        Map<String,Integer> fields = createFieldsIndexMap(lines.get(0), TSV_DELIM);

        for(String line : lines.subList(1, lines.size()))
        {
            String[] values = line.split(TSV_DELIM, -1);
            sigEtiologies.put(values[fields.get("signature")], values[fields.get("etiology")]);
        }

        return sigEtiologies;
    }

    public static void populateBucketMap(final Map<String,Integer> bucketNameIndexMap)
    {
        populateBucketMap(bucketNameIndexMap, null);
    }

    public static void populateBucketMap(final Map<String,Integer> bucketNameIndexMap, @Nullable final List<String> bucketNames)
    {
        char[] refBases = {'C', 'T'};
        char[] bases = {'A','C', 'G', 'T'};
        int index = 0;

        for(int i = 0; i < refBases.length; ++i)
        {
            char ref = refBases[i];

            for(int j = 0; j < bases.length; ++j)
            {
                char alt = bases[j];

                if(ref != alt)
                {
                    String baseChange = String.format("%c>%c", ref, alt);

                    for (int k = 0; k < bases.length; ++k)
                    {
                        char before = bases[k];

                        for (int l = 0; l < bases.length; ++l)
                        {
                            char after = bases[l];

                            String context = String.format("%c%c%c", before, ref, after);

                            String bucketName = baseChange + "_" + context;

                            if(bucketNames != null)
                                bucketNames.add(bucketName);

                            bucketNameIndexMap.put(bucketName, index);
                            ++index;
                        }
                    }
                }
            }
        }
    }

    public static String contextFromVariant(final SmallVariant variant)
    {
        return variantContext(variant.ref(), variant.alt(), variant.trinucleotideContext());
    }

    public static String variantContext(final String ref, final String alt, final String trinucContext)
    {
        // convert base change to standard set and the context accordingly
        if(ref.charAt(0) == 'A' || ref.charAt(0) == 'G')
        {
            return String.format("%c>%c_%c%c%c",
                    swapDnaBase(ref.charAt(0)), swapDnaBase(alt.charAt(0)),
                    swapDnaBase(trinucContext.charAt(2)), swapDnaBase(trinucContext.charAt(1)), swapDnaBase(trinucContext.charAt(0)));
        }
        else
        {
            return String.format("%c>%c_%s",ref.charAt(0), alt.charAt(0), trinucContext);
        }
    }
}
