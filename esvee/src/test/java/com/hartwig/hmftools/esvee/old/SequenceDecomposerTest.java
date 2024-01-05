package com.hartwig.hmftools.esvee.old;

public class SequenceDecomposerTest
{
    /* CHASHA FIXME
    @Test
    public void canDecompose()
    {
        final String bases = "GCCTGGCTATATATATATATATTTTTTTTTTTTTTTTTTTTAGTAGAGATGGGGTTTCACTGTGTTG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(4);
        assertTrue(result.toString()).isEqualTo("[GCCTGGC, TAx7, Tx20, AGTAGAGATGGGGTTTCACTGTGTTG]");
    }

    @Test
    public void canDecomposeRepeatStart()
    {
        final String bases = "TTTTTTTTTTTTTTTTTTTTTATATATATATATAAGTAGAGATGGGGTTTCACTGTGTTG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[Tx21, ATx6, AAGTAGAGATGGGGTTTCACTGTGTTG]");
    }

    @Test
    public void canDecomposeRepeat2Start()
    {
        final String bases = "TATATATATATATATTTTTTTTTTTTTTTTTTTTAGTAGAGATGGGGTTTCACTGTGTTG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[TAx7, Tx20, AGTAGAGATGGGGTTTCACTGTGTTG]");
    }

    @Test
    public void canDecomposeRepeat3Start()
    {
        final String bases = "TATTATTATTATTATTATTATTTTTTTTTTTTTTTTTTTTTAGTAGAGATGGGGTTTCACTGTGTTG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[TATx7, Tx20, AGTAGAGATGGGGTTTCACTGTGTTG]");
    }

    @Test
    public void canDecomposeRepeat3WithDouble()
    {
        final String bases = "TTATTATTATTATTATTATTATTTTTTTTTTTTTTTTTTTTTAGTAGAGATGGGGTTTCACTGTGTTG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[TTAx7, Tx21, AGTAGAGATGGGGTTTCACTGTGTTG]");
    }

    @Test
    public void canDecomposeRepeat3Middle()
    {
        final String bases = "ATACGGCTAATCGTATTATTATTATTATTATTATTATATCGGCTAATCG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[ATACGGCTAATCG, TATx8, ATCGGCTAATCG]");
    }

    @Test
    public void canDecomposeRepeat4Start()
    {
        final String bases = "TATCTATCTATCTATCTATCTATCTATCCCCCCCCCCCCCCCCCCCCCAGTAGAGATGGGGTTTCACTGTGTTG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[TATCx7, Cx20, AGTAGAGATGGGGTTTCACTGTGTTG]");
    }

    @Test
    public void canDecomposeRepeat4AndRepair()
    {
        final String bases = "TATCTATCTATCTATCTATCTATCTATCTATACCCCCCCCCCCCCCCCCCCCAGTAGAGATGGGGTTTCACTGTGTTG";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(4);
        assertTrue(result.toString()).isEqualTo("[TATCx7, TATA, Cx20, AGTAGAGATGGGGTTTCACTGTGTTG]");
    }

    @Test
    public void canDecomposeRepeat4Whole()
    {
        final String bases = "TATCTATCTATCTATCTATCTATCTATC";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(1);
        assertTrue(result.toString()).isEqualTo("[TATCx7]");
    }

    @Test
    public void canDecomposeRepeat4Partial()
    {
        final String bases = "TATCTATCTATCTATCTATCTATCTAT";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(2);
        assertTrue(result.toString()).isEqualTo("[TATCx6, TAT]");
    }

    @Ignore
    @Test
    public void canDecomposeRepeat6With5BaseRepeatFront()
    {
        final String bases = "TATCAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATTATC";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[TATC, AAAAATx9, TATC]");
    }

    @Ignore
    @Test
    public void canDecomposeRepeat6With5BaseRepeatBack()
    {
        final String bases = "TATCTAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATATC";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[TATC, TAAAAAx9, TATC]");
    }

    @Test
    public void canDecomposeRepeat6With5BaseRepeatWhole()
    {
        final String bases = "AAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAAT";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(1);
        assertTrue(result.toString()).isEqualTo("[AAAAATx9]");
    }

    @Ignore
    @Test
    public void canDecomposeRepeat6With5BaseRepeatEdge()
    {
        final String bases = "GCAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAGC";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(3);
        assertTrue(result.toString()).isEqualTo("[GC, AAATAAx9, GC]");
    }

    @Test
    public void doesNotConfuseRepeat1WithRepeat6()
    {
        final String bases = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);
        final List<SequenceDecomposer.Node> result = SequenceDecomposer.decompose(bases.getBytes(), quals);
        assertTrue(result).hasSize(1);
        assertTrue(result.toString()).isEqualTo("[Ax54]");
    }
    */
}