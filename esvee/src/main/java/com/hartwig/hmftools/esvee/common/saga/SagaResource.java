package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SagaResource
{
    private final String mFastaPath;
    private final List<SagaAssembly> mAssemblies;
    private final SAMSequenceDictionary mSequenceDictionary;

    public SagaResource(final String fastaPath)
    {
        mFastaPath = fastaPath;
        mAssemblies = new ArrayList<>(100_000);
        mSequenceDictionary = loadAssemblies(fastaPath, mAssemblies);
    }

    public String bwaIndexImagePath()
    {
        return mFastaPath + ".img";
    }

    public SAMSequenceDictionary samDict()
    {
        return mSequenceDictionary;
    }

    public List<SagaAssembly> assemblies()
    {
        return mAssemblies;
    }

    private static SAMSequenceDictionary loadAssemblies(final String fastaPath, List<SagaAssembly> assemblies)
    {
        SAMSequenceDictionary samDict;
        try(IndexedFastaSequenceFile fasta = loadFasta(fastaPath))
        {
            samDict = fasta.getSequenceDictionary();
            samDict.getSequences().stream()
                    .map(seq -> SagaAssembly.fromFastaLabel(seq.getSequenceName(), seq.getSequenceLength()))
                    .forEach(assemblies::add);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to load SAGA resource", e);
        }

        // Variants are keyed by ID later, so ensure the IDs are unique.
        long uniqueVariantIdCount = assemblies.stream().map(assembly -> assembly.variant().id()).distinct().count();
        if(uniqueVariantIdCount != assemblies.size())
        {
            throw new RuntimeException("Duplicate variant IDs in SAGA resource");
        }

        SV_LOGGER.debug("loaded {} SAGA variant assemblies from file({})", assemblies.size(), fastaPath);

        return samDict;
    }

    private static IndexedFastaSequenceFile loadFasta(final String fastaPath) throws IOException
    {
        return new IndexedFastaSequenceFile(new File(fastaPath));
    }
}
