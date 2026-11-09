package com.hartwig.hmftools.amber.purity;

import java.util.ArrayList;
import java.util.List;

public class LocalMaximaFinder<T extends Score>
{
    private enum Direction
    {
        UP,
        DOWN,
        FLAT
    }

    private static int FIRST_PEAK_ITEM_COUNT = 3;
    private static double FIRST_PEAK_MIN_RATIO = 1.25;
    private static double PEAK_VS_TROUGH_MIN_RATIO = 0.03;

    private final List<T> mScores;

    public LocalMaximaFinder(final List<T> scores)
    {
        mScores = scores;
    }

    public List<T> maxima()
    {
        List<T> maxima = new ArrayList<>();
        T firstNonZero = mScores.stream().filter(score -> score.score() != 0.0).findFirst().orElse(null);
        if(firstNonZero == null)
        {
            return maxima;
        }

        Direction previousDirection = null;
        Direction currentDirection = null;
        T previousScore = null;

        double firstNonZeroScore = 0;
        double lastTroughScore = -1;
        int nonZeroScoreCount = 0;

        for(T score : mScores)
        {
            if(previousScore != null)
            {
                if(score.score() > previousScore.score())
                {
                    currentDirection = Direction.UP;
                }
                else if(score.score() < previousScore.score())
                {
                    currentDirection = Direction.DOWN;
                }
                else
                {
                    currentDirection = Direction.FLAT;
                }

                if(currentDirection == Direction.DOWN)
                {
                    // down after peak
                    boolean addPeak = false;

                    if(previousDirection == Direction.UP)
                    {
                        // avoid creating a peak close to another at the start
                        if(nonZeroScoreCount <= FIRST_PEAK_ITEM_COUNT)
                        {
                            if(previousScore.score() >= firstNonZeroScore * FIRST_PEAK_MIN_RATIO)
                                addPeak = true;
                        }
                        else
                        {
                            addPeak = true;
                        }
                    }

                    if(addPeak)
                    {
                        if(lastTroughScore > 0 && !maxima.isEmpty())
                        {
                            int lastIndex = maxima.size() - 1;
                            double lastPeak = maxima.get(lastIndex).score();

                            double peakToTroughRatio = (lastPeak - lastTroughScore) / lastTroughScore;

                            if(peakToTroughRatio < PEAK_VS_TROUGH_MIN_RATIO)
                                maxima.remove(lastIndex);
                        }

                        maxima.add(previousScore);

                        lastTroughScore = -1;
                    }
                }
                else if(previousDirection == Direction.DOWN && (currentDirection == Direction.UP || currentDirection == Direction.FLAT))
                {
                    lastTroughScore = previousScore.score();
                }
            }

            if(score.score() > 0)
            {
                ++nonZeroScoreCount;

                if(firstNonZeroScore == 0)
                    firstNonZeroScore = score.score();
            }

            if(currentDirection != Direction.FLAT)
            {
                previousDirection = currentDirection;
                previousScore = score;
            }
        }

        maxima.remove(firstNonZero);
        return maxima;
    }
}
