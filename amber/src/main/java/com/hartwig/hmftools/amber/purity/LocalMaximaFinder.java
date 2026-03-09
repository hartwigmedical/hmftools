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

    private final List<T> Scores;

    public LocalMaximaFinder(final List<T> scores)
    {
        Scores = scores;
    }

    public List<T> maxima()
    {
        List<T> maxima = new ArrayList<>();
        T firstNonZero = Scores.stream().filter(score -> score.score() != 0.0).findFirst().orElse(null);
        if(firstNonZero == null)
        {
            return maxima;
        }
        Direction previousDirection = null;
        Direction currentDirection = null;
        T previousScore = null;
        for(T score : Scores)
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
                        maxima.add(previousScore);
                    }
                }
            }

            previousDirection = currentDirection;
            previousScore = score;
        }
        maxima.remove(firstNonZero);
        return maxima;
    }
}
