# The example function below keeps track of the opponent's history and plays whatever the opponent played two plays ago. It is not a very good player so you will need to change the code to pass the challenge.

def player(prev_play, opponent_history=[]):
    opponent_history.append(prev_play)
    guess = "S"

    most_frequent = max(set(opponent_history), key=opponent_history.count)

    if most_frequent == "R": guess = "P"
    if most_frequent == "P": guess = "S"
    if most_frequent == "S": guess = "R"
    last_ten = opponent_history[-10:]
    if most_frequent == '':
        most_frequent = "S"

    ideal_response = {'P': 'S', 'R': 'P', 'S': 'R'}
    guess = ideal_response[most_frequent]
    return guess
