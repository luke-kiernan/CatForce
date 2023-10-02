#include "LifeAPI.h"
static const bool CP_DEBUG = false;
static const int MAX_CP_CHAIN_LEN = -1;

// build up state from a single cell
// via next state in chain = unmoved | (moved shifted by (xShift, yShift))

// record chain using 4n+2 integers:
// n entries that record the states: unmovedIndex, movedIndex, xShift, yShift
// one last entry that records: xPostShift, yPostShift

std::vector<int>  ComputeCPChain(const LifeState & desired){ //, bool debug){
  std::vector<LifeState> chainedStates;
  std::vector<LifeState> shiftsInsideDesired;
  std::vector<int> cpChain;
  LifeState oneCell;
  oneCell.Set(0,0);
  chainedStates.push_back(oneCell);
  shiftsInsideDesired.push_back(desired);

  // weird edge cases
  assert(!desired.IsEmpty());
  if (desired.GetPop() == 1){
    auto onCell = desired.FirstOn();
    cpChain.push_back(onCell.first);
    cpChain.push_back(onCell.second);
    return cpChain;
  }

  // greedy algorithm. choose the next state in chain to maximize
  // area, tie-breaking by number of shifts that fall inside desired region.
  bool done = false;
  unsigned count = 0;
  while(!done){
    ++count;
    unsigned bestArea = 0;
    unsigned numShiftsOfBestInDesired = 0;
    LifeState best;
    for(unsigned i = 0; i < 4; ++i)
      cpChain.push_back(0); // write dummy data to end of chain (will overwrite with best)
    for(int i = chainedStates.size()-1; i >= 0; --i){
      if (chainedStates[i].GetPop() +(chainedStates.end()-1)->GetPop() < bestArea || done)
        break;
      for(int j = chainedStates.size()-1; j >= i; --j){
        if (chainedStates[i].GetPop() + chainedStates[j].GetPop() < bestArea || done)
          break;
        // relShifts = [shifts of moved in desired] - [shifts of unmoved in desired]
        // i is unmoved, j is moved
        LifeState shiftsOfUnmoved = shiftsInsideDesired[i];
        shiftsOfUnmoved.Transform(Rotate180OddBoth);
        LifeState relShifts = shiftsInsideDesired[j].Convolve(shiftsOfUnmoved);
        while(!relShifts.IsEmpty()){
          auto onCell = relShifts.FirstOn();
          relShifts.Erase(onCell.first, onCell.second);
          LifeState unioned = chainedStates[i].UnionWithShifted(chainedStates[j],
                    onCell.first, onCell.second);
          unsigned numShiftsInDesired = 0;
          if(unioned.GetPop() == bestArea)
            numShiftsInDesired = desired.MatchLive(unioned).GetPop();
          // did we find a better one?
          if(unioned.GetPop() > bestArea || 
                  numShiftsInDesired > numShiftsOfBestInDesired){
            bestArea = unioned.GetPop();
            best = unioned;
            const std::array<int, 4> data({i,j, onCell.first, onCell.second});
            std::copy(data.begin(), data.end(), cpChain.end()-4);
            if(bestArea == desired.GetPop()){
              done = true;
              break;
            } else if (numShiftsInDesired == 0)
              numShiftsOfBestInDesired = desired.MatchLive(unioned).GetPop();
            else
              numShiftsOfBestInDesired = numShiftsInDesired;
          }
        }
      }
    }
    if(CP_DEBUG){
      std::cout << "best for state " << count << std::endl;
      best.Print();
    }
    chainedStates.push_back(best);
    shiftsInsideDesired.emplace_back(desired.MatchLive(best));
    if(done){
      auto postShift = (shiftsInsideDesired.end()-1)->FirstOn();
      cpChain.push_back(postShift.first);
      cpChain.push_back(postShift.second);
    }
    //assert(count <= desired.GetPop());
  }
  if ( MAX_CP_CHAIN_LEN > 0 && cpChain.size() > MAX_CP_CHAIN_LEN){
    std::cout << "max copy paste chain length too small";
    std::cout << ": increase to at least " << cpChain.size() << std::endl;
    assert(false);
  }
  return cpChain;
}

LifeState ApplyCPChainTo(const LifeState & inState,
                        const std::vector<int> & chain,
                        bool unionType){ //, bool debug){
  std::vector<LifeState> chainedStates;
  chainedStates.reserve((chain.size()-2)/4+1);
  chainedStates.push_back(inState);
  if (CP_DEBUG){
    std::cout << "state 0: input state" << std::endl;
    inState.Print();
  }
  for(unsigned i = 1; i <= (chain.size()-2)/4; ++i){
    unsigned startBlock = 4*(i-1);
    const LifeState & unmoved = chainedStates[chain[startBlock]];
    const LifeState & moved = chainedStates[chain[startBlock+1]];
    if(unionType)
      chainedStates.emplace_back(unmoved.UnionWithShifted(moved,
            chain[startBlock+2], chain[startBlock+3]));
    else // [reversing of direction built into IntersectWithShifted]
      chainedStates.emplace_back(unmoved.IntersectWithShifted(moved,
            chain[startBlock+2], chain[startBlock+3]));
    if (CP_DEBUG){
      std::string verb = unionType ? std::string("union") : std::string("intersect");
      std::cout << "state " << i << ": state " << chain[startBlock];
      std::cout << " " << verb << " state " << chain[startBlock+1] << " shifted by ";
      std::cout << chain[startBlock+2] << " " << chain[startBlock+3] << std::endl;
      (chainedStates.end()-1)->Print();
    }
  }
  if(unionType)
    chainedStates[chainedStates.size()-1].Move(*(chain.end()-2),*(chain.end()-1));
  else // need to reverse the direction.
    chainedStates[chainedStates.size()-1].Move(64-*(chain.end()-2),64-*(chain.end()-1));
  if (CP_DEBUG)
    std::cout << "shifting final state by (" << *(chain.end()-2) << ", " << *(chain.end()-1) << ")." << std::endl;
  return chainedStates[chainedStates.size()-1];
}