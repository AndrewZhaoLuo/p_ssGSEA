import numpy.random as random
import timeit

'''
Simulates shuffling a deck of cards, checking whether a certain sample of cards is on "top" of the deck.
The "top" of the deck is considered to be the index 0.

num_cards           =   number of cards to include in the deck, they will be numbered from [0, num_cards]. should be
                        positive

cutoff_index        =   cards being sampled should, after shuffling in the deck, be above or equal to this index.
                        should be non-negative and less than num_cards

sample_size         =   the number of cards to make sure are over the cutoff_index. this will be done by selecting
                        and tracking the cards numbered from [0, sample_size) when shuffled. should be less than or
                        equal to num_cards

num_trials          =   the number of shuffles and trackings to run. should be positive

returns a tuple (num_success, num_fails) where num_success is the number of shuffles where all tracked cards were
over the cutoff index and num_fails is the opposite
'''
def shufflesim(num_cards, cutoff_index, sample_size, num_trials):
    num_success = 0
    num_fails = 0

    for trial in range(0, num_trials):
        shuffledDeck = random.permutation(num_cards)
        check = lambda deck: sum([1 for x in range(0, cutoff_index + 1) if deck[x] < sample_size]) == sample_size

        #check to see if cards aren't where its supposed to be (below cutoff)
        if check(shuffledDeck):
            num_success += 1
        else:
            num_fails += 1

    return (num_success, num_fails)

if __name__ == "__main__":
    #enviroment variables
    big_tab = "\t\t\t\t"
    small_tab = "\t"

    cards = 52
    cutoff_index = int((cards-1) / 2)
    sample_size = 10

    trials = 100000

    #start the shuffling!
    start = timeit.default_timer()

    num_success, num_fails = shufflesim(cards, cutoff_index, sample_size, trials)
    print("+:")
    print(small_tab + str(num_success))
    print("-:")
    print(small_tab + str(num_fails))
    print()

    end = timeit.default_timer()

    print(small_tab + " Cards:")
    print(big_tab + str(cards))
    print(small_tab + " Cutoff:")
    print(big_tab + str(cutoff_index))
    print(small_tab + " Samples:")
    print(big_tab + str(sample_size))
    print(small_tab + " Trials:")
    print(big_tab + str(trials))
    print()

    print("Took:" + str(end-start) + "s")


