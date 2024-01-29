package com.hartwig.hmftools.bamtools.bamtofastq.collection;

import java.util.Collection;
import java.util.concurrent.atomic.AtomicReference;

import org.jetbrains.annotations.Nullable;
import org.pcollections.ConsPStack;
import org.pcollections.PStack;

// TODO NEXT: TEST
// TODO: Move to HMF commons.
public class AtomicPStack<T>
{
    private final AtomicReference<PStack<T>> mStackRef;

    public AtomicPStack()
    {
        mStackRef = new AtomicReference<>(ConsPStack.empty());
    }

    public AtomicPStack(final Collection<T> elements)
    {
        PStack<T> stack = ConsPStack.empty();
        for(T element : elements)
        {
            stack = stack.plus(element);
        }

        mStackRef = new AtomicReference<>(stack);
    }

    public void push(T element)
    {
        mStackRef.getAndUpdate(stack -> stack.plus(element));
    }

    @Nullable
    public T pop()
    {
        PStack<T> currentStack = mStackRef.getAndUpdate(stack ->
        {
            if(stack.isEmpty())
            {
                return stack;
            }

            return stack.minus(0);
        });

        if(currentStack.isEmpty())
        {
            return null;
        }

        return currentStack.get(0);
    }

    @Nullable
    public T peek()
    {
        PStack<T> currentStack = mStackRef.get();
        if(currentStack.isEmpty())
        {
            return null;
        }

        return currentStack.get(0);
    }

    public int size()
    {
        return mStackRef.get().size();
    }

    public boolean isEmpty()
    {
        return mStackRef.get().isEmpty();
    }
}
