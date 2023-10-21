#include "LifeAPI.h"
#include <string.h>
#include <cassert>
static const bool CP_DEBUG = false;
static const int MAX_CP_CHAIN_LEN = -1;
// could switch to arrays of int8_t's.

// build up state from a single cell
// via next state in chain = unmoved | (moved shifted by (xShift, yShift))

// record chain using 4n+2 integers:
// n entries that record the states: unmovedIndex, movedIndex, xShift, yShift
// one last entry that records: xPostShift, yPostShift

enum OperationType {Union, Intersection};

std::vector<int> ComputeCPChain(const LifeState & desired){
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
    for (unsigned m = 0; m < 4; ++m)
      cpChain.push_back(0); // write dummy data to end, will overwrite later.
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
    if ( MAX_CP_CHAIN_LEN > 0 && cpChain.size() > MAX_CP_CHAIN_LEN){
      std::cout << "max copy paste chain length too small";
      std::cout << ": increase to at least " << cpChain.size() << std::endl;
      assert(false);
    }
    if(done){ // not sure if this is necessary. absolute position shouldn't matter.
    // plus maybe this should be the bounding box, not the first on.
      auto postShift = (shiftsInsideDesired.end()-1)->FirstOn();
      cpChain.push_back(postShift.first);
      cpChain.push_back(postShift.second);
    }
  }
  return cpChain;
}

LifeState ApplyCPChainTo(const LifeState & inState,
                        const std::vector<int> & chain,
                        OperationType type){
                        // unionType is true for convolution, false for pattern matching
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
    if(type == OperationType::Union)
      chainedStates.emplace_back(unmoved.UnionWithShifted(moved,
            chain[startBlock+2], chain[startBlock+3]));
    else // [reversing of direction built into IntersectWithShifted]
      chainedStates.emplace_back(unmoved.IntersectWithShifted(moved,
            chain[startBlock+2], chain[startBlock+3]));
    if (CP_DEBUG){
      std::string verb = (type == OperationType::Union) 
                      ? std::string("union") : std::string("intersect");
      std::cout << "state " << i << ": state " << int(chain[startBlock]);
      std::cout << " " << verb << " state " << int(chain[startBlock+1]) << " shifted by ";
      std::cout << int(chain[startBlock+2]) << " " << int(chain[startBlock+3]) << std::endl;
      (chainedStates.end()-1)->Print();
    }
  }
  unsigned chain_len = chain.size();
  if(type == OperationType::Union)
    chainedStates[chainedStates.size()-1].Move(chain[chain_len-2],chain[chain_len-1]);
  else // need to reverse the direction.
    chainedStates[chainedStates.size()-1].Move(64-chain[chain_len-2],64-chain[chain_len-1]);
  if (CP_DEBUG){
    std::cout << "shifting final state by (" << int(chain[chain_len-2]);
    std::cout << ", " << int(chain[chain_len-1]) << ")." << std::endl;
  }
  return chainedStates[chainedStates.size()-1];
}

class CP_Target{
  private:
    std::vector<int> wanted_chain;
    std::vector<int> unwanted_chain;
  public:
    CP_Target(const LifeState & desired){
      wanted_chain = ComputeCPChain(desired);
      LifeState halo = desired.ZOI() & ~desired;
      unwanted_chain = ComputeCPChain(halo);
    }
    CP_Target() = default;
    // compatible with match-survive parameter
    // (Useful for conduit searching: eliminates results where output
    // object runs straight into a catalyst.)
    bool MatchLiveAndDead(const LifeState & workspace, const LifeState & catalysts,
                          int maxJunk, int matchSurvive){
      unsigned int chain_len = wanted_chain.size();
      if (chain_len == 0)
        return true;
      std::array<SymmetryTransform, 8> transforms = {Identity, Rotate90, Rotate180OddBoth,
                      Rotate270, ReflectAcrossX, ReflectAcrossY, ReflectAcrossYeqX,
                      ReflectAcrossYeqNegXP1};
      LifeState matchIn = workspace & ~catalysts;
      for(auto transf : transforms){
        LifeState temp_pos = matchIn;
        LifeState temp_neg = ~matchIn;
        temp_pos.Transform(transf);
        temp_neg.Transform(transf);
        LifeState matches  = ApplyCPChainTo(temp_pos, wanted_chain, OperationType::Intersection) &
                              ApplyCPChainTo(temp_neg, unwanted_chain, OperationType::Intersection);
        if (!matches.IsEmpty()){
          LifeState junk;
          if (matchSurvive > 0){
            LifeState matchedAdvanced = ApplyCPChainTo(matches, wanted_chain, OperationType::Union);
            matchedAdvanced.Step(matchSurvive);
            temp_pos |= catalysts;
            temp_pos.Step(matchSurvive);
            if (!(temp_pos.Contains(matchedAdvanced)))
              continue;
            junk = temp_pos & ~matchedAdvanced & ~catalysts;
          } else{
            junk = temp_pos & ~(ApplyCPChainTo(matches, wanted_chain, OperationType::Union));
          }
          if (maxJunk < 0 || junk.GetPop() < maxJunk)
            return true;
        }
      }
      return false;
    }
    void Print(){
      for (unsigned i = 0; i < wanted_chain.size(); ++i){
        std::cout << wanted_chain[i] << ( (i % 4 == 3) ? "; " : " ");
        if ((i+1) % 20 == 0) std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    
};