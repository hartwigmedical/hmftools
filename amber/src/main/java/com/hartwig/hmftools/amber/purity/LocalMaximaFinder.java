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

                if(previousDirection == Direction.UP || previousDirection == Direction.FLAT)
                {
                    if(currentDirection == Direction.DOWN)
                    {
                        // avoid creating a peak close to another at the start
                        if(nonZeroScoreCount <= 3)
                        {
                            if(previousScore.score() >= firstNonZeroScore * 1.25)
                                maxima.add(previousScore);
                        }
                        else
                        {
                            maxima.add(previousScore);
                        }
                    }
                }
            }

            if(score.score() > 0)
            {
                ++nonZeroScoreCount;

                if(firstNonZeroScore == 0)
                    firstNonZeroScore = score.score();
            }

            previousDirection = currentDirection;
            previousScore = score;
        }

        maxima.remove(firstNonZero);
        return maxima;
    }
}
