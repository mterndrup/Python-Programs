# The example function below keeps track of the opponent's history and plays whatever the opponent played two plays ago. It is not a very good player so you will need to change the code to pass the challenge.

def player(prev_play, opponent_history=[]):
    opponent_history.append(prev_play)
    guess = "R"
    guesses = "RPS"

    most_frequent = max(set(opponent_history), key=opponent_history.count)

    if most_frequent == "R": guess = "P"
    if most_frequent == "P": guess = "S"
    if most_frequent == "S": guess = "R"

    return guess
