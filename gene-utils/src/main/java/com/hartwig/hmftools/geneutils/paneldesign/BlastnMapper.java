package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.paneldesign.BlastnResult.INVALID_RESULT;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.BLASTN_WORD_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.MIN_BLAST_ALIGNMENT_LENGTH;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.blastn.BlastnMatch;
import com.hartwig.hmftools.common.blastn.BlastnRunner;

public class BlastnMapper
{
    private final PanelConfig mConfig;

    private final BlastnRunner mBlastnRunner;
    private final BlastnCache mResultsCache;

    public BlastnMapper(final PanelConfig config)
    {
        mConfig = config;

        if(config.BlastTool != null && config.BlastDb != null)
        {
            mBlastnRunner = new BlastnRunner.Builder()
                    .withTask("megablast")
                    .withPrefix(config.OutputPrefix)
                    .withBlastDir(config.BlastTool)
                    .withBlastDb(config.BlastDb)
                    .withOutputDir(config.OutputDir)
                    .withKeepOutput(true)
                    .withWordSize(BLASTN_WORD_SIZE)
                    .withSubjectBestHit(true)
                    .withNumThreads(config.Threads)
                    .build();
        }
        else
        {
            mBlastnRunner = null;
        }

        mResultsCache = new BlastnCache(mConfig.BlastCacheFile);
    }

    public Multimap<Integer,BlastnMatch> mapSequences(final Map<Integer,String> querySequences)
    {
        if(mBlastnRunner == null)
        {
            // TODO: simulate a passing result if skipping BlastN calls
            return ArrayListMultimap.create();
        }

        return mBlastnRunner.run(querySequences);
    }

    public List<BlastnResult> mapSequences(final List<String> sequences)
    {
        List<BlastnResult> results = Lists.newArrayListWithCapacity(sequences.size());

        if(mBlastnRunner == null)
        {
            if(!mConfig.SkipBlast)
                return Collections.emptyList();

            for(int i = 0; i < sequences.size(); ++i)
            {
                String sequence = sequences.get(i);
                results.add(new BlastnResult(sequence, 1, 1));
            }

            return results;
        }

        // initialise results which actual or cached results may then override
        for(int i = 0; i < sequences.size(); ++i)
        {
            results.add(INVALID_RESULT);
        }

        Map<Integer,String> sequencesMap = Maps.newHashMap();

        for(int i = 0; i < sequences.size(); ++i)
        {
            String sequence = sequences.get(i);

            if(mResultsCache.enabled())
            {
                BlastnResult result = mResultsCache.findResult(sequence);

                if(result != null)
                {
                    results.set(i, result);
                    continue;
                }
            }

            sequencesMap.put(i, sequence);
        }

        if(sequencesMap.isEmpty())
            return results;

        Multimap<Integer,BlastnMatch> sequenceResults = mBlastnRunner.run(sequencesMap);

        List<BlastnResult> newResults = Lists.newArrayList();

        for(Integer i : sequencesMap.keySet())
        {
            String sequence = sequences.get(i);

            Collection<BlastnMatch> matches = sequenceResults.get(i);

            if(matches.isEmpty())
            {
                GU_LOGGER.debug("sequence({}) returned no blastn results", sequence);
                BlastnResult result = new BlastnResult(sequence, 0, 0);
                results.set(i, result);
                newResults.add(result);
                continue;
            }

            double sumBitScore = 0;

            // process all the matches and sum up the bit score, but remove the one with best match
            Optional<BlastnMatch> bestMatch = matches.stream()
                    .filter(x -> isPrimaryBlastnMatch(x))
                    .max(Comparator.comparing(BlastnMatch::getBitScore));

            if(bestMatch.isPresent())
            {
                sumBitScore -= bestMatch.get().getBitScore();
            }

            for(BlastnMatch m : matches)
            {
                if(m.getAlignmentLength() >= MIN_BLAST_ALIGNMENT_LENGTH && isPrimaryBlastnMatch(m))
                {
                    sumBitScore += m.getBitScore();
                }
            }

            BlastnResult result = new BlastnResult(sequence, sumBitScore, matches.size());
            results.set(i, result);
            newResults.add(result);
        }

        persistResult(newResults);

        return results;
    }

    private void persistResult(final List<BlastnResult> results)
    {
        if(mConfig.BlastCacheFile == null)
            return;

        results.forEach(x -> mResultsCache.addCacheEntry(x));
    }

    public void onComplete()
    {
        mResultsCache.writeCache();
    }

    public static boolean isPrimaryBlastnMatch(BlastnMatch match)
    {
        return match.getSubjectTitle().contains("Primary Assembly") ||
                match.getSubjectTitle().contains("unlocalized genomic scaffold") ||
                match.getSubjectTitle().contains("unplaced genomic scaffold");
    }
}
