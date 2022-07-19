package com.hartwig.hmftools.cdr3.layout

import com.hartwig.hmftools.cdr3.ReadKey
import htsjdk.samtools.util.SequenceUtil
import org.apache.logging.log4j.LogManager
import java.util.IdentityHashMap

class LayoutTree(val minBaseQuality: Byte, val minOverlapLength: Int)
{
    class Read (
        // have a source that allows us to refer back to where this comes from
        val source: Any,
        val readKey: ReadKey,
        val sequence: String,
        val baseQualities: ByteArray,
        val layoutPosition: Int) // position relative to layout start

    // each node has a base and their associated count
    class Node(val level: Int, var parent: Node?)
    {
        internal var base: Char = UNKNOWN_BASE
        internal var highQualBase: Char = UNKNOWN_BASE
        internal var highQualityCount: Int = 0
        val reads: MutableList<Read> = ArrayList()
        val children: MutableList<Node> = ArrayList()
    }

    // root does not correspond to any node level
    val root: Node = Node(-1, null)
    private val levelNodes: MutableList<MutableList<Node>> = ArrayList()

    val numLevels: Int get() { return levelNodes.size }

    fun matchesNode(node: Node, base: Char, baseQuality: Byte) : Boolean
    {
        return node.base == base || baseQuality < minBaseQuality || node.highQualBase == UNKNOWN_BASE || base == UNKNOWN_BASE
    }

    fun getOrCreateLevelNodes(level: Int) : MutableList<Node>
    {
        while (levelNodes.size <= level)
        {
            levelNodes.add(ArrayList())
        }
        return levelNodes[level]
    }

    fun tryAddRead(read: Read) : Boolean
    {
        if (read.layoutPosition == 0)
        {
            // we can start at the root
            addReadToBranch(read, root)
            return true
        }
        else
        {
            val level = read.layoutPosition
            if (levelNodes.size <= level + minOverlapLength)
            {
                // we are not up to this level yet, cannot possibly build from here
                return false
            }

            // we start from all the potential parents
            var added = false
            val branches: List<Node> = levelNodes[level - 1]

            for (branch in branches)
            {
                if (overlapsBranch(read, branch))
                {
                    addReadToBranch(read, branch)
                    added = true
                }
            }
            return added
        }
    }

    // we want to test that this read overlaps this branch
    private fun overlapsBranch(read: Read, branchRoot: Node) : Boolean
    {
        assert(read.layoutPosition > 0)
        assert(branchRoot.level == read.layoutPosition - 1)

        var currentNode = branchRoot
        for (i in 0 until Math.min(minOverlapLength, read.sequence.length))
        {
            val baseQual: Byte = read.baseQualities[i]
            val base: Char = read.sequence[i]
            val highQualBase: Char = if (baseQual >= minBaseQuality) base else UNKNOWN_BASE
            val childNode: Node? = chooseChildNode(currentNode, base, highQualBase)

            if (childNode == null)
                return false

            currentNode = childNode
        }
        return true
    }

    private fun addReadToBranch(read: Read, branchRoot: Node)
    {
        val leafNodes = ArrayList<Node>()

        findMatchingLeafNodes(read, branchRoot, -1, leafNodes)

        if (leafNodes.isNotEmpty())
        {
            val leaf = leafNodes[0]

            assert(leaf.children.isEmpty())

            var current = leaf

            // add all necessary child nodes
            for (i in leaf.level - read.layoutPosition + 1 until read.sequence.length)
            {
                val baseQual: Byte = read.baseQualities[i]
                val base: Char = read.sequence[i]
                val highQualBase: Char = if (baseQual >= minBaseQuality) base else UNKNOWN_BASE

                current = createChild(current)
                current.base = base

                if (highQualBase != UNKNOWN_BASE)
                {
                    current.highQualBase = highQualBase
                    current.highQualityCount++
                }
            }

            current.reads.add(read)

            // update all nodes from leaf to root
            current = leaf
            while (current != root)
            {
                val i = read.layoutPosition + current.level
                val baseQual: Byte = read.baseQualities[i]
                val base: Char = read.sequence[i]
                val highQualBase: Char = if (baseQual >= minBaseQuality) base else UNKNOWN_BASE

                if (highQualBase != UNKNOWN_BASE)
                {
                    current.highQualityCount++
                    if (current.highQualBase == UNKNOWN_BASE)
                    {
                        current.base = highQualBase
                        current.highQualBase = highQualBase

                        // now we have to try to merge it with parents, which is also a bit of work
                        tryMergeWithSibling(current)
                    }
                }
                current = current.parent!!
            }

            return
        }

        // add a new branch
        var currentNode = branchRoot
        for (i in read.sequence.indices)
        {
            val level = read.layoutPosition + i
            assert(currentNode.level + 1 == level)
            val baseQual: Byte = read.baseQualities[i]
            val base: Char = read.sequence[i]
            val highQualBase: Char = if (baseQual >= minBaseQuality) base else UNKNOWN_BASE

            var child: Node? = null

            if (highQualBase == UNKNOWN_BASE)
            {
                // this is the secret sauce here. We always create a new child node if
                // we have a low qual base. This allow subsequence reads to match here even
                // if this base does not match
                currentNode = createChild(currentNode)
                currentNode.base = base
                continue
            }

            for (c in currentNode.children)
            {
                if (c.highQualBase == highQualBase)
                {
                    // we found a node we can assign to
                    c.highQualityCount++
                    child = c
                    break
                }
            }

            if (child == null)
            {
                // create a new one
                child = createChild(currentNode)
                child.base = base
                child.highQualBase = highQualBase
                child.highQualityCount++
            }

            currentNode = child
        }

        // this is end
        currentNode.reads.add(read)
    }

    private fun findMatchingLeafNodes(read: Read, node: Node, i: Int, leafNodes: MutableList<Node>)
    {
        if (i >= 0)
        {
            val level = read.layoutPosition + i
            assert(node.level == level)
            val baseQual: Byte = read.baseQualities[i]
            val base: Char = read.sequence[i]

            if (!matchesNode(node, base, baseQual))
            {
                return
            }

            if (node.children.isEmpty())
            {
                leafNodes.add(node)
                return
            }
        }

        for (child in node.children)
        {
            findMatchingLeafNodes(read, child, i + 1, leafNodes)
        }
    }

    // we have a node that we just went from unknown / low qual base to high qual base
    // we try to merge it with sibling
    private fun tryMergeWithSibling(node: Node)
    {
        if (node.parent == null)
            return

        for (sibling in node.parent!!.children)
        {
            if (sibling === node)
                continue

            if (sibling.highQualBase == node.highQualBase)
            {
                // we can merge them
                mergeNodes(sibling, node)
                break
            }
        }
    }

    private fun mergeNodes(nodeTo: Node, nodeFrom: Node)
    {
        nodeTo.highQualityCount += nodeFrom.highQualityCount
        nodeTo.parent!!.children.removeIf({ n -> n === nodeFrom })

        for (r in nodeFrom.reads)
        {
            if (!nodeTo.reads.any({ o -> o === r }))
            {
                nodeTo.reads.add(r)
            }
        }

        for (childFrom in nodeFrom.children)
        {
            var merged = false
            for (childTo in nodeTo.children)
            {
                if (childFrom.highQualBase != UNKNOWN_BASE &&
                    childFrom.highQualBase == childTo.highQualBase)
                {
                    mergeNodes(childTo, childFrom)
                    merged = true
                    break
                }
            }
            if (!merged)
            {
                nodeTo.children.add(childFrom)
                childFrom.parent = nodeTo
                continue
            }
        }
    }

    fun chooseChildNode(parent: Node, base: Char, highQualBase: Char) : Node?
    {
        for (node in parent.children)
        {
            if (node.base == base)
            {
                // we found a node we can assign to
                return node
            }
        }

        if (highQualBase == UNKNOWN_BASE)
        {
            // can't really match anything, just match first one. Should not be too much of issue
            if (parent.children.size > 0)
            {
                return parent.children[0]
            }
            // if there is no children then can't do anything
            return null
        }
        else
        {
            // if we cannot find one then see if any is unknown base match
            for (node in parent.children)
            {
                if (node.highQualBase == UNKNOWN_BASE)
                {
                    // we found a node we can assign to
                    return node
                }
            }
            return null
        }
    }

    fun createChild(parent: Node) : Node
    {
        val lvl = parent.level + 1
        val child = Node(parent.level + 1, parent)
        parent.children.add(child)
        getOrCreateLevelNodes(lvl).add(child)
        return child
    }

    // trim away branches that have low support
    fun trimBranches(minHighQualSupport: Int)
    {
        // start from root
        val toVisit = ArrayList<Node>()

        for (node: Node in toVisit)
        {
            // can't merge branch if only one child
            if (node.children.size <= 1)
            {
                continue
            }

            for (child: Node in node.children)
            {
                // if any of the child has low support, try to
                // put it into the other child
                if (child.highQualityCount < minHighQualSupport)
                {
                    // we go through the siblings and see if we can reassign this
                    // read to the other ones
                    for (siblings: Node in node.children)
                    {
                        if (child == siblings)
                            continue

                    }
                }
            }
        }
    }

    fun buildReadLayouts(useReverseComp: Boolean = false) : List<ReadLayout>
    {
        val layouts = ArrayList<ReadLayout>()

        // depth first traversal, we visit every leaf node and build a layout from it
        val leafNodes = ArrayList<Node>()
        val levelChildIndex = IntArray(levelNodes.size, { 0 })
        var current: Node? = root

        while (current != null)
        {
            if (current.children.isEmpty())
            {
                // this is a leaf
                leafNodes.add(current)
                // go back up one level
                current = current.parent
            }
            else
            {
                // which child are we up to?
                val childIndex = levelChildIndex[current.level + 1]
                if (childIndex == current.children.size)
                {
                    // finished with this node
                    levelChildIndex[current.level + 1] = 0
                    current = current.parent
                }
                else
                {
                    levelChildIndex[current.level + 1] += 1
                    current = current.children[childIndex]
                }
            }
        }

        sLogger.info("layout tree: found {} leaf nodes", leafNodes.size)

        // we want to avoid creating too many read objects. We want to share them
        val readToLayoutReadMap = IdentityHashMap<Read, ReadLayout.Read>()

        for (leaf in leafNodes)
        {
            layouts.add(buildReadLayout(leaf, readToLayoutReadMap, useReverseComp))
        }

        sLogger.info("built {} layouts from layout tree", layouts.size)

        return layouts
    }

    // create a layout from the leaf node
    fun buildReadLayout(leafNode: Node, readToLayoutReadMap: IdentityHashMap<Read, ReadLayout.Read>, useReverseComp: Boolean = false) : ReadLayout
    {
        sLogger.debug("building layout from node")

        if (leafNode.children.isNotEmpty())
        {
            throw IllegalArgumentException("cannot build layout from non leaf node")
        }

        val reads = ArrayList<ReadLayout.Read>()
        val readKeys = HashSet<ReadKey>()

        val layoutLength = leafNode.level + 1

        var current: Node? = leafNode
        while (current != null && current != root)
        {
            // building it from leaf, we get all the reads
            for (r in current.reads)
            {
                // NOTE: it would only make sense to add read that ends here

                if (readKeys.add(r.readKey))
                {
                    var layoutRead = readToLayoutReadMap[r]

                    if (layoutRead == null)
                    {
                        // we need to create layout read from this read
                        if (useReverseComp)
                        {
                            layoutRead = ReadLayout.Read(r.source,
                                r.readKey, SequenceUtil.reverseComplement(r.sequence),
                                r.baseQualities.reversed().toByteArray(),
                                r.sequence.length - r.layoutPosition - 1)
                        }
                        else
                        {
                            layoutRead = ReadLayout.Read(r.source,
                                r.readKey, r.sequence,
                                r.baseQualities,
                                r.layoutPosition)
                        }
                        readToLayoutReadMap[r] = layoutRead
                    }
                    reads.add(layoutRead)
                }
            }
            current = current.parent
        }

        val layout = ReadLayout()

        for (r in reads)
        {
            layout.addRead(r, minBaseQuality)
        }

        sLogger.debug("layout: num reads({}) seq({}) support({})", layout.reads.size, layout.consensusSequence(), layout.highQualSupportString())

        return layout
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(LayoutTree::class.java)
        const val UNKNOWN_BASE: Char = 'N'
    }
}