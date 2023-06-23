package com.hartwig.hmftools.cider.layout

import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList

//
// This class creates a tree to help layout the reads from the most probable to least
//
// The problem it is trying to solve is
//
// Say we have 3 reads that support the following:
//   ACTG
// next read we get is
//   ACTGG
// but the last G is actually an error, such that we get
//   ACTGG   <--- error
//   ACTGT
//   ACTGT
//   ACTGT
//
// Using the old method of sequentially building layouts, we cannot recover from this G error,
// and all the reads up to G will be assigned to this layout. Which gives incorrect impression that
// this G is highly supported.
// We need to be able to recover from the error and assign those reads to the ACTGT layout instead
//
// To solve this problem, we build a tree, which looks like this when we have the first 3 ACTG reads:
//
// root
//  |
//  A3-C3-T3-G3
//
// the number denotes the read count at each node.
// When we get the ACTGG read, the G is added as a child:
// root
//  |
//  A4-C4-T4-G4-G1
//
// And when the 3 ACTGT reads are added, the last T do not match, so it is added as a branch node:
//
// root
//  |
//  A7-C7-T7-G7-G1
//            \
//             T3   <---- the last T do not match the existing G
//
// Note that there is a constraint where there must be reads that ends at the G at the 4th position.
// We cannot create a branch if we have the following reads: ATACTGG, ACTGT. It may look like the last
// base can be either G or T, but since ATACTGG is a full read, the only possibility is ATACTGG, and
// ATACTGT is not possible, since it will involve taking first 6 bases of first read and 1 base of 2nd read.
//
// When we extract the layout from the tree, we start from the leaf nodes that follows the path of
// highest support. In this case, the leaf node at T3 will be built first since it has highest support,
// and all reads except the one with error will be assigned to this layout.
//
// There is another problem we have to consider, from the above tree:
//
// root
//  |
//  A7-C7-T7-G7-G1
//            \
//             T3
//
// Say we build it out further, and they both follow same path, and we want to add 3 reads of CCA, i.e.
//
// root
//  |
//  A7-C7-T7-G7-G2-C1-C1-A1
//             \
//              T3-C3-C1-A1
//
// When a new read of CCA is added, we want to put it at the branch with highest support, i.e. the bottom
// T branch instead. Even though CCA matches with both.
//
// To solve this problem, we keep track of all the nodes that start on a given level. And when a new read
// is added that can match multiple branches, we use the one with highest support
//
// Now consider another read of TCCT, we can add it to the above tree and get:
//
// root
//  |
//  A7-C7-T7-G7-G2-C1-C1-A1
//             \
//              T4-C4-C2-A1
//                     \
//                      T1
// However, if we have another read of ACCA, we cannot add it to the tree. Since it does not overlap with
// any of the existing branches. We have to create a totally new branch from the root where the
// earlier nodes are empty, i.e.
//
// root
//  | \
//  | A7-C7-T7-G7-G2-C1-C1-A1
//  |            \
//  |             T4-C4-C2-A1
//  |                     \
//  |                      T1
//  |-------------A1-C1-C1-A1
// Next consider we have a read CCA, where should it be added? It should be added to the top branch
// since top branch has higher support. We have to test both branches when it can go to more than one place
//
// Another important note is that we should not be adding branches to nodes that are already highly supported
// Once a node has a certain number of base support they should be sealed.
//
class LayoutTree(val minBaseQuality: Byte, val minOverlapLength: Int, var minSupportToSealNode: Int)
{
    class Read (
        // have a source that allows us to refer back to where this comes from
        val source: Any,
        val sequence: String,
        val baseQualities: ByteArray,
        val alignedPosition: Int) // where in the read is the aligned position

    // each node has a base and their associated count
    class Node(var position: Int, var parent: Node?) : Comparable<Node>
    {
        internal var base: Char = UNKNOWN_BASE
        internal var highQualBase: Char = UNKNOWN_BASE
        internal var count: Int = 0
        internal var highQualityCount: Int = 0
        val reads: MutableList<Read> = ArrayList()
        val children: MutableList<Node> = ArrayList()

        // support is the high quality count + a little bit of the non
        // high quality count.
        val support : Double get()
        {
            require(highQualityCount <= count)
            return highQualityCount + (count - highQualityCount) * LOW_QUAL_BASE_FACTOR
        }

        override fun compareTo(other: Node): Int
        {
            return compare(this, other)
        }

        companion object
        {
            // compare the support, if they have the same support then compare the parents
            private fun compare(node1: Node, node2: Node): Int
            {
                var n1: Node? = node1
                var n2: Node? = node2

                while (n1 != null && n2 != null)
                {
                    if (n1 === n2)
                        return 0
                    if (n1.support < n2.support)
                        return -1
                    if (n1.support > n2.support)
                        return 1
                    n1 = n1.parent
                    n2 = n2.parent
                }
                return 0
            }
        }
    }

    // root does not correspond to any read position
    val root: Node = Node(-1, null)
    var alignedPosition: Int = Int.MIN_VALUE

    // we store the nodes of each level. This is used primarily to allow us to
    // jump to a level where a read starts.
    private val levelNodes: MutableList<MutableList<Node>> = ArrayList()

    val numLevels: Int get() { return levelNodes.size }

    fun toLayoutPosition(readAlignedPos: Int) : Int
    {
        return alignedPosition - readAlignedPos
    }

    fun matchesNode(node: Node, base: Char, baseQuality: Byte) : Boolean
    {
        return node.base == base || baseQuality < minBaseQuality || node.highQualBase == UNKNOWN_BASE || base == UNKNOWN_BASE
    }

    fun createChild(parent: Node, lvl: Int) : Node
    {
        require(parent === root || lvl == parent.position + 1)

        // check that we are not adding children to a sealed node
        require(parent === root || !isNodeSealed(parent))

        val child = Node(lvl, parent)
        parent.children.add(child)
        getOrCreateLevelNodes(lvl).add(child)
        return child
    }

    fun getOrCreateLevelNodes(level: Int) : MutableList<Node>
    {
        while (levelNodes.size <= level)
        {
            levelNodes.add(ArrayList())
        }
        return levelNodes[level]
    }

    fun addToNode(node: Node, base: Char, baseQuality: Byte)
    {
        val highQualBase: Char = if (baseQuality >= minBaseQuality) base else UNKNOWN_BASE

        if (highQualBase != UNKNOWN_BASE)
        {
            require(node.highQualBase == UNKNOWN_BASE || node.highQualBase == highQualBase)
            node.highQualityCount++
            if (node.highQualBase == UNKNOWN_BASE)
            {
                node.base = highQualBase
                node.highQualBase = highQualBase
            }
        }
        if (base != UNKNOWN_BASE && node.base == base)
        {
            node.count++
        }
    }

    fun isNodeSealed(node: Node) : Boolean
    {
        return node.children.any { n -> n.highQualityCount >= minSupportToSealNode }
    }

    private fun checkUpdateAlignedPosition(readAlignedPos: Int)
    {
        if (readAlignedPos <= alignedPosition)
        {
            return
        }

        val diff = readAlignedPos - alignedPosition

        // insert levels in the level nodes list
        if (levelNodes.isNotEmpty())
        {
            for (nodeList in levelNodes)
            {
                nodeList.forEach({ n -> n.position += diff })
            }

            for (i in 0 until diff)
            {
                levelNodes.add(0, ArrayList())
            }

            sLogger.warn("WARNING: adding reads in wrong order. Reads should be added from largest alignedPosition to smallest")
        }

        alignedPosition = readAlignedPos
    }

    fun tryAddRead(read: Read) : Boolean
    {
        checkUpdateAlignedPosition(read.alignedPosition)

        val readLayoutPos = toLayoutPosition(read.alignedPosition)

        if (levelNodes.size >= readLayoutPos + minOverlapLength)
        {
            // find the branch that overlaps with this read, and choose
            // the one with highest support
            var bestBranch: Node? = null

            for (branch in levelNodes[readLayoutPos])
            {
                if ((bestBranch == null || branch > bestBranch) && overlapsBranch(read, branch))
                {
                    // this branch has higher support and overlaps with the read
                    bestBranch = branch
                }
            }

            if (bestBranch != null)
            {
                require(bestBranch.parent!!.children.contains(bestBranch))
                require(levelNodes[bestBranch.position].contains(bestBranch))
                addReadToBranch(read, bestBranch)
                return true
            }
        }

        // add this read to root
        // create a new one
        val branchRoot = Node(readLayoutPos, root)
        root.children.add(branchRoot)
        getOrCreateLevelNodes(readLayoutPos).add(branchRoot)
        branchRoot.base = read.sequence[0]
        branchRoot.highQualBase = if (read.baseQualities[0] >= minBaseQuality) branchRoot.base else UNKNOWN_BASE
        addReadToBranch(read, branchRoot)
        return true
    }

    // branchRoot is the place where the first base of the read can
    // be added.
    // We test to see if the first minOverlapLength bases matches with this branch
    // The definition of matching is either the base is low quality, or the bases must be the
    // same if both are high quality.
    // Once the first min overlap length matches, we test the rest of the sequence to determine
    // if it requires branching any nodes that are sealed
    private fun overlapsBranch(read: Read, branchRoot: Node) : Boolean
    {
        val readLayoutPos = toLayoutPosition(read.alignedPosition)
        require(branchRoot !== root)
        require(readLayoutPos >= 0)
        require(branchRoot.position == readLayoutPos)

        var currentNode = branchRoot

        if (!matchesNode(branchRoot, read.sequence[0], read.baseQualities[0]))
        {
            return false
        }

        // in the following we check for 2 things. Firstly if the min overlap has been satisfied.
        // And if minoverlap has been satisfied, we want to check if the rest of the sequences
        // can be added to the tree without branching a node that is sealed
        loop@ for (i in 1 until read.sequence.length)
        {
            val baseQual: Byte = read.baseQualities[i]
            val base: Char = read.sequence[i]
            val highQualBase: Char = if (baseQual >= minBaseQuality) base else UNKNOWN_BASE

            if (i >= minOverlapLength && currentNode.children.isEmpty())
            {
                // we reached the end of this tree, need to extend it, that is totally fine
                break
            }

            for (node in currentNode.children)
            {
                if (node.base == base)
                {
                    // we found a node we can assign to
                    currentNode = node
                    continue@loop
                }
            }

            if (highQualBase == UNKNOWN_BASE)
            {
                // we should have already protected against this case with previous check
                assert(i < minOverlapLength || currentNode.children.isNotEmpty())

                // can't really match anything, match the one with highest number of support
                // if there is no node then return false
                currentNode = currentNode.children.maxByOrNull { n -> n.support } ?: return false
                continue@loop
            }
            // if we cannot find one then see if any is unknown base match
            for (node in currentNode.children)
            {
                if (node.highQualBase == UNKNOWN_BASE)
                {
                    // we found a node we can assign to
                    currentNode = node
                    continue@loop
                }
            }

            if (i >= minOverlapLength)
            {
                // if we get to here we need to see if we can add a branch off this node. We can
                // only add branch if the node is not sealed
                return !isNodeSealed(currentNode)
            }

            // cannot find a child
            return false
        }

        return true
    }

    private fun addReadToBranch(read: Read, branchRoot: Node)
    {
        // we have already established that this reads overlaps with this branch
        val readLayoutPos = toLayoutPosition(read.alignedPosition)

        // follow this branch until a mismatch occur and we create a new subbranch from there.
        // node that this logic has to match the overlap function
        var currentNode = branchRoot
        for (i in 1 until read.sequence.length)
        {
            val level = readLayoutPos + i
            assert(currentNode.position + 1 == level)
            val baseQual: Byte = read.baseQualities[i]
            val base: Char = read.sequence[i]
            val highQualBase: Char = if (baseQual >= minBaseQuality) base else UNKNOWN_BASE

            var child: Node? = null

            for (c in currentNode.children)
            {
                if (c.base == base)
                {
                    // we found a node we can assign to
                    child = c
                    addToNode(child, base, baseQual)
                    break
                }
            }

            if (child == null)
            {
                if (highQualBase == UNKNOWN_BASE)
                {
                    // if it is unknown, we just put it into the branch with highest support
                    // this is not entirely correct, but seems to work just fine
                    child = currentNode.children.maxByOrNull { n -> n.support }
                }
                else
                {
                    // we find a branch with unknown base
                    child = currentNode.children.find { n -> n.highQualBase == UNKNOWN_BASE }
                }
            }

            if (child == null)
            {
                // create a new one
                child = createChild(currentNode, level)
                child.base = base
                addToNode(child, base, baseQual)
            }

            currentNode = child
        }

        // this is end
        currentNode.reads.add(read)
    }

    private fun removeNode(node: Node)
    {
        // remove it from the nodes level list
        assert(levelNodes[node.position].contains(node))
        val removed = levelNodes[node.position].remove(node)

        require(removed)
    }

    // visit all the nodes in depth first manner
    fun depthFirstVisit(startNode: Node, visitor: (Node) -> Unit)
    {
        // depth first traversal
        val indexStack = Stack<Int>()
        var current: Node = startNode
        indexStack.push(0)

        while (true)
        {
            val childIndex: Int = indexStack.peek()
            if (childIndex == 0)
            {
                // first time we visit this node
                // NOTE: we must allow visitor to modify the children at this
                //       point. Never assume children list stays the same
                visitor(current)
            }

            if (childIndex == current.children.size)
            {
                // finished with this node
                if (current == startNode)
                {
                    // we are back to the start node
                    return
                }
                indexStack.pop()
                current = current.parent!!
                continue
            }
            // update current node index
            indexStack.push(indexStack.pop() + 1)

            // add child node index
            indexStack.push(0)
            current = current.children[childIndex]
        }
    }

    // for every node we sort the children descending order by highQualityCount
    fun sortNodeChildren()
    {
        sLogger.debug("sorting children")

        depthFirstVisit(root, { n -> n.children.sortByDescending( { c -> c.support } )})
        sLogger.debug("finished sorting children")
    }

    // We try to reassign reads by removing branches that are not the top branch and try to re add the
    // reads to the tree. This will progressively seal more and more nodes until we cannot reassign
    // any reads
    fun reassignReads()
    {
        val original = minSupportToSealNode

        for (i in 0 until 10) // try maximum of 10 rounds
        {
            val readsToReassign = ArrayList<Read>()

            // depth first traversal, we visit every leaf node and build a layout from it
            // NOTE this is not quite correct as we remove the reads but the counts are not removed
            // from the upstream nodes
            depthFirstVisit(root) { node: Node ->
                if (node !== root && node.children.size > 1 && isNodeSealed(node))
                {
                    // we choose just 1 child and remove the rest
                    val childToKeep: Node = node.children.maxByOrNull { n -> n.support }!!

                    for (c in node.children)
                    {
                        if (c !== childToKeep)
                        {
                            // put all reads into the reassign list, also remove all
                            depthFirstVisit(c) { n: Node ->
                                readsToReassign.addAll(n.reads)
                                removeNode(n)
                            }
                        }
                    }
                    node.children.clear()
                    node.children.add(childToKeep)
                }
            }

            if (readsToReassign.isNotEmpty())
            {
                sLogger.debug("reassigning {} reads", readsToReassign.size)

                // now re add all those reads
                for (r in readsToReassign)
                {
                    tryAddRead(r)
                }
            }
            else if (minSupportToSealNode <= 1)
            {
                break
            }

            minSupportToSealNode = (minSupportToSealNode - 1).coerceAtLeast(1)
        }

        minSupportToSealNode = original
    }

    // we do depth first traversal, by following the path with highest support
    // after a read is used we remove it. This means that the reads will tend
    // to go towards the highest supported branch.
    fun buildReadLayouts(layoutReadCreate: (read: Read) -> ReadLayout.Read) : List<ReadLayout>
    {
        // first have to sort every node children from highest to lowest support
        sortNodeChildren()
        reassignReads()

        val layouts = ArrayList<ReadLayout>()

        val leafNodes = ArrayList<Node>()

        // depth first traversal, we visit every leaf node and build a layout from it
        depthFirstVisit(root) { node: Node ->
            if (node !== root && node.children.isEmpty())
            {
                // this is a leaf node
                leafNodes.add(node)
            }
        }

        sLogger.trace("layout tree: found {} leaf nodes", leafNodes.size)

        // a set to keep track of reads that we have already used
        val usedReadSet: MutableSet<Read> = Collections.newSetFromMap(IdentityHashMap())

        for (leaf in leafNodes)
        {
            layouts.add(buildLayoutFromLeaf(leaf, layoutReadCreate, usedReadSet))
        }

        sLogger.trace("built {} layouts from layout tree", layouts.size)

        return layouts
    }

    // create a layout from the leaf node
    fun buildLayoutFromLeaf(leafNode: Node,
                        layoutReadCreate: (read: Read) -> ReadLayout.Read,
                        usedReadSet: MutableSet<Read>) : ReadLayout
    {
        sLogger.trace("building layout from node")

        if (leafNode.children.isNotEmpty())
        {
            throw IllegalArgumentException("cannot build layout from non leaf node")
        }

        val reads = ArrayList<ReadLayout.Read>()

        var current: Node? = leafNode
        while (current != null && current != root)
        {
            // building it from leaf, we get all the reads
            for (r in current.reads)
            {
                // NOTE one way or the other we must ensure the read is not added
                // multiple times into a layout. Since a read exists in many levels we
                // need to have a check
                if (!usedReadSet.add(r))
                {
                    continue
                }

                // we need to create layout read from this read
                val layoutRead: ReadLayout.Read = layoutReadCreate(r)
                reads.add(layoutRead)
            }
            current = current.parent
        }

        // we want to build the reads the reverse way, i.e. from leftmost
        // to rightmost
        reads.reverse()

        val layout = ReadLayout()

        for (r in reads)
        {
            layout.addRead(r, minBaseQuality)
        }

        sLogger.trace("layout: num reads({}) seq({}) support({})", layout.reads.size, layout.consensusSequence(), layout.highQualSupportString())

        return layout
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(LayoutTree::class.java)
        const val UNKNOWN_BASE: Char = 'N'
        const val LOW_QUAL_BASE_FACTOR: Double = 0.1
    }
}
