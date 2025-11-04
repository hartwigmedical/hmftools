package com.hartwig.hmftools.common.genome.refgenome;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public class CachedRefGenome implements RefGenomeInterface
{
    private static final int DEFAULT_BLOCK_SIZE = 1_000;
    private static final int DEFAULT_MAX_CACHED_BLOCKS = 3;

    private static class Block
    {
        public final String Chromosome;
        public final int Idx;

        public Block(final String chromosome, int idx)
        {
            Chromosome = chromosome;
            Idx = idx;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
                return true;

            if(!(o instanceof final Block block))
                return false;

            return Idx == block.Idx && Objects.equals(Chromosome, block.Chromosome);
        }

        @Override
        public int hashCode()
        {
            return Chromosome.hashCode() + 31 * Idx;
        }
    }

    private class BlockCache extends LinkedHashMap<Block, byte[]>
    {
        public BlockCache()
        {
            super(mMaxCachedBlocks + 1, 0.75f, true);
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry<Block, byte[]> eldest)
        {
            return size() > mMaxCachedBlocks;
        }
    }

    private final RefGenomeInterface mRefGenome;
    private final ThreadLocal<BlockCache> mThreadBlockCaches = new ThreadLocal<BlockCache>()
    {
        @Override protected BlockCache initialValue() { return new BlockCache(); }
    };
    private final int mBlockSize;
    private final int mMaxCachedBlocks;

    public CachedRefGenome(final RefGenomeInterface refGenome, int blockSize, int maxCachedBlocks)
    {
        mRefGenome = refGenome;
        mBlockSize = blockSize;
        mMaxCachedBlocks = maxCachedBlocks;
    }

    public CachedRefGenome(final RefGenomeInterface refGenome)
    {
        this(refGenome, DEFAULT_BLOCK_SIZE, DEFAULT_MAX_CACHED_BLOCKS);
    }

    public final RefGenomeInterface refGenome() { return mRefGenome; }

    @Override
    public String getBaseString(final String chromosome, final int posStart, final int posEnd)
    {
        return new String(getBases(chromosome, posStart, posEnd));
    }

    @Override
    public String getBaseString(final String chromosome, final List<int[]> baseRanges)
    {
        StringBuilder bases = new StringBuilder();
        for(int[] range : baseRanges)
        {
            int start = range[0];
            int end = range[1];
            bases.append(getBaseString(chromosome, start, end));
        }

        return bases.toString();
    }

    @Override
    public int getChromosomeLength(final String chromosome)
    {
        return mRefGenome.getChromosomeLength(chromosome);
    }

    @Override
    public byte[] getBases(final String chromosome, final int posStart, final int posEnd)
    {
        int chrStartPos = oneBasedIndexing() ? 1 : 0;
        int chrEndPos = getChromosomeLength(chromosome) - 1 + chrStartPos;
        if(posStart < chrStartPos || posEnd > chrEndPos)
            throw new IllegalArgumentException(format("Requested ref genome region out of bounds: %s:%d-%d", chromosome, posStart, posEnd));

        byte[] bases = new byte[posEnd - posStart + 1];
        int basesIdx = 0;
        int startBlockIdx = (posStart - chrStartPos) / mBlockSize;
        int endBlockIdx = (posEnd - chrStartPos) / mBlockSize;
        for(int i = startBlockIdx; i <= endBlockIdx; i++)
        {
            byte[] block = getBlock(chromosome, i);
            int blockStart = i * mBlockSize + chrStartPos;
            int blockEnd = min(blockStart + mBlockSize - 1, chrEndPos);
            int startIdx = max(posStart, blockStart) - blockStart;
            int endIdx = min(posEnd, blockEnd) - blockStart;
            int basesCopied = endIdx - startIdx + 1;
            arraycopy(block, startIdx, bases, basesIdx, basesCopied);
            basesIdx += basesCopied;
        }

        return bases;
    }

    @Override
    public byte getBase(final String chromosome, int pos)
    {
        int chrStartPos = oneBasedIndexing() ? 1 : 0;
        int chrEndPos = getChromosomeLength(chromosome) - 1 + chrStartPos;
        if(pos < chrStartPos || pos > chrEndPos)
            throw new IllegalArgumentException(format("Requested ref genome base out of bounds: %s:%d", chromosome, pos));

        int blockIdx = (pos - chrStartPos) / mBlockSize;
        byte[] block = getBlock(chromosome, blockIdx);
        int blockStart = blockIdx * mBlockSize + chrStartPos;
        int idx = pos - blockStart;
        return block[idx];
    }

    @Override
    public Map<String, Integer> chromosomeLengths()
    {
        return mRefGenome.chromosomeLengths();
    }

    private byte[] getBlock(final String chromosome, int blockIdx)
    {
        BlockCache blockCache = mThreadBlockCaches.get();
        Block block = new Block(chromosome, blockIdx);
        byte[] bases = blockCache.getOrDefault(block, null);
        if(bases != null)
            return bases;

        int chrStartPos = oneBasedIndexing() ? 1 : 0;
        int chrEndPos = getChromosomeLength(chromosome) - 1 + chrStartPos;
        int blockStart = blockIdx * mBlockSize + chrStartPos;
        int blockEnd = min(blockStart + mBlockSize - 1, chrEndPos);
        bases = mRefGenome.getBases(chromosome, blockStart, blockEnd);
        blockCache.put(block, bases);
        return bases;
    }

    @Override
    public boolean oneBasedIndexing() { return mRefGenome.oneBasedIndexing(); }
}
